#include "comp_llpmatrix.h"
//#include <mex.h>
#include <map>

struct LessdPoint{ 
	public: LessdPoint(){};
	bool operator()(const dPoint& p1, const dPoint& p2){ 
		for(int i = 0; i < p1.dimension(); i ++){
			if(p1[i] < p2[i]){ return true; }
			else if((p1[i] > p2[i])){ return false;}
		}
		return false; 
	} 
};



void generate_laplace_matrix_sparse_matlab_dim2(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS)
{
//	fprintf(stdout, "generate_laplace_matrix_sparse_matlab_dim2\n");

	typedef CGAL::Search_traits_d<KCd> Traits;
	typedef CGAL::Fair<Traits> Fair;
	typedef CGAL::Kd_tree<Traits, Fair> Tree;
	typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
	typedef KCd::FT FT;

	const unsigned int tdim = 2;

	//cout<<"h: "<<h<<" rho: "<<rho<<endl;
	double hh = h * h;
	double nmfactor = M_PI * hh * hh / 4;
	//cout<<"nmfactor: "<<nmfactor<<endl;
//	fprintf(stdout, "nmfactor: %f\n", nmfactor);

	unsigned int np = pcloud.p_count();
	vector<dPoint> points;
	map<const dPoint, int, LessdPoint> pt2index;
	for(unsigned int i = 0; i < np; i ++){
		points.push_back( pcloud.point(i).coord() );
		pt2index.insert( make_pair(pcloud.point(i).coord(), i) );
	}
	Fair fair(10);
	Tree tree(points.begin(), points.end(), fair);
	points.clear();

	
	//cout<<"\n Come Here GLMSMD2 1.0\n";


	double eps = 1e-3 * h * rho;
	vector<dPoint> neighbor_pts;
	vector<int>	neighbor_indices;
	vector<dVector> tspace;
	Triangulation_2 tr;


	//cout<<"\n Come Here GLMSMD2 2.0\n";

	for(unsigned int i = 0; i < np; i ++){
		
		//cout<<"i: "<<i<<"\r";
//		fprintf(stdout, "i: %d\r", i);

		dPoint pt = pcloud.point(i).coord();
		
		//search the points within distance h * rho;
		neighbor_pts.clear();
		Fuzzy_sphere fs(pt, h * rho, FT(eps));
		tree.search(back_insert_iterator<vector<dPoint> >(neighbor_pts),  fs);

		//look up the indices of the neighboring points
		neighbor_indices.clear();

		bool find = false;
	   for(vector<dPoint>::iterator iter = neighbor_pts.begin(); iter != neighbor_pts.end(); iter ++){
  			//cout<<*iter<<endl;
			map<const dPoint, int, LessdPoint>::iterator iter_pt2index = pt2index.find(*iter);
			if( iter_pt2index != pt2index.end() ){
				neighbor_indices.push_back(iter_pt2index->second);
				if(iter_pt2index->second == (int)i){
					find = true;
				}
			}
			else{
				//cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
				fprintf(stderr,"Fail to find the index of the neighboring points\n");
			}
   	}
		assert(neighbor_indices.size() >= 3 && find);

		//cout<<i<<"th neighbors: ";
		//for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
		//	cout<< neighbor_indices[j]<<" ";
		//}
		//cout<<endl;

		
		//estimate normal
		tspace.clear();
		est_tangent_space_lapack(pt, neighbor_pts, h, tdim, tspace);
		
		//cout<<"pt: "<<pt<<endl;
		//cout<<"tspace: "<<tspace[0]<<", "<<tspace[1]<<endl;

		tr.clear();	
		for(unsigned int j = 0; j < neighbor_pts.size(); j ++){
			dVector vec = neighbor_pts[j] - pt;

			//cout<<vec<<endl;
			//cout<<neighbor_indices[j]<<endl;
			double x = CGAL::to_double(vec * tspace[0]);
			double y = CGAL::to_double(vec * tspace[1]);
			Vertex_handle_2d vh = tr.insert( Point_2d(x, y) );
			assert(vh != NULL);
			if(vh->get_origin_id() >= 0){
				//cout<<"warning: containing duplicated points"<<endl;
				fprintf(stderr,"Containing duplicated points\n");

			}
			vh->set_origin_id( neighbor_indices[j] );
		}

		map<unsigned int, double > vareas;
		for(Triangulation_2::Vertex_iterator viter = tr.vertices_begin(); viter != tr.vertices_end(); viter ++){
			vareas[viter->get_origin_id()] = 0;
		}
		for(Triangulation_2::Face_iterator fiter = tr.faces_begin(); fiter != tr.faces_end(); fiter ++){
			Triangle_2d triangle( fiter->vertex(0)->point(), fiter->vertex(1)->point(), fiter->vertex(2)->point() );
			double farea = fabs( CGAL::to_double(triangle.area()) );
			vareas[fiter->vertex(0)->get_origin_id()] += farea / 3;
			vareas[fiter->vertex(1)->get_origin_id()] += farea / 3;
			vareas[fiter->vertex(2)->get_origin_id()] += farea / 3;
		}

		//--------------------------------------------------
		// //debug
		// cout<<"after computing area: "<<endl;
		// for(map<unsigned int, double>::iterator iter = vareas.begin(); iter != vareas.end(); iter ++){
		// 	cout<< "vid: "<<iter->first << " area: "<< iter->second <<endl; 
		// }
		//-------------------------------------------------- 
		
		double totalweight = 0;

		for(Triangulation_2::Vertex_iterator viter = tr.vertices_begin(); viter != tr.vertices_end(); viter ++){
			unsigned int vid = viter->get_origin_id();
			if( vid == i ){
				assert( (viter->point() - Point_2d(0, 0)) * (viter->point() - Point_2d(0, 0)) < FT(1e-5) );
				continue;
			}
			//double weight = exp(-CGAL::to_double( (viter->point() - Point_2d(0, 0)) * (viter->point() - Point_2d(0, 0))) / hh) * ( 4.0 / (M_PI * hh * hh) );
			double sqdist = CGAL::to_double( (viter->point() - Point_2d(0, 0)) * (viter->point() - Point_2d(0, 0)) );
			double weight = exp(-sqdist / hh) / nmfactor;

			//cout<<"vid: "<<vid<<" dist: "<< sqdist <<" weight: "<<weight<<endl;

			II.push_back(i + 1);
			JJ.push_back(vid + 1);
			SS.push_back(weight * vareas[vid]);
			//IIfout<<i + 1<<" ";
			//JJfout<<vid + 1<<" ";
			//SSfout<<weight * vareas[vid]<<" ";
			totalweight -= weight * vareas[vid];
		}
		II.push_back(i + 1);
		JJ.push_back(i + 1);
		SS.push_back(totalweight);
		//IIfout<<i + 1<<" ";
		//JJfout<<i + 1<<" ";
		//SSfout<<totalweight<<" ";
	}
	fprintf(stdout, "\n");
}

void generate_laplace_matrix_sparse_matlab_dim3(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS)
{
	//cout<<"generate_laplace_matrix_sparse_matlab_dim3"<<endl;
	fprintf(stdout, "generate_laplace_matrix_sparse_matlab_dim3\n");

	typedef CGAL::Search_traits_d<KCd> Traits;
	typedef CGAL::Fair<Traits> Fair;
	typedef CGAL::Kd_tree<Traits, Fair> Tree;
	typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
	typedef KCd::FT FT;


	const unsigned int tdim = 3;

	//cout<<"h: "<<h<<" rho: "<<rho<<endl;
	double hh = h * h;
	double nmfactor = sqrt(M_PI * hh) * M_PI * hh * hh / 4;
	//cout<<"nmfactor: "<<nmfactor<<endl;
	fprintf(stdout, "nmfactor: %f\n", nmfactor);

	unsigned int np = pcloud.p_count();
	vector<dPoint> points;
	map<const dPoint, int, LessdPoint> pt2index;
	for(unsigned int i = 0; i < np; i ++){
		points.push_back( pcloud.point(i).coord() );
		pt2index.insert( make_pair(pcloud.point(i).coord(), i) );
	}
	Fair fair(10);
	Tree tree(points.begin(), points.end(), fair);
	points.clear();

	double eps = 1e-3 * h * rho;
	vector<dPoint> neighbor_pts;
	vector<int>	neighbor_indices;
	vector<dVector> tspace;
	Triangulation tr;

	for(unsigned int i = 0; i < np; i ++){

		//cout<<"i: "<<i<<endl;
		fprintf(stdout, "i: %d\r", i);
		
		dPoint pt = pcloud.point(i).coord();
		
		//search the points within distance h * rho;
		neighbor_pts.clear();
		Fuzzy_sphere fs(pt, h * rho, FT(eps));
		tree.search(back_insert_iterator<vector<dPoint> >(neighbor_pts),  fs);

		//look up the indices of the neighboring points
		neighbor_indices.clear();

		bool find = false;
	   for(vector<dPoint>::iterator iter = neighbor_pts.begin(); iter != neighbor_pts.end(); iter ++){
  			//cout<<*iter<<endl;
			map<const dPoint, int, LessdPoint>::iterator iter_pt2index = pt2index.find(*iter);
			if( iter_pt2index != pt2index.end() ){
				neighbor_indices.push_back(iter_pt2index->second);
				if(iter_pt2index->second == (int)i){
					find = true;
				}
			}
			else{
				//cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
				fprintf(stderr,"Fail to find the index of the neighboring points\n");
			}
   	}
		assert(neighbor_indices.size() >= 3 && find);

		//cout<<i<<"th neighbors: ";
		//for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
		//	cout<< neighbor_indices[j]<<" ";
		//}
		//cout<<endl;

		
		//estimate normal
		tspace.clear();
		est_tangent_space_lapack(pt, neighbor_pts, h, tdim, tspace);
		
		//pcloud.vertex(i).set_normal(normal);	
		//cout<<"pt: "<<pt<<endl;
		//cout<<"tspace: "<<tspace[0]<<", "<<tspace[1]<<endl;
		//getchar();

		//project the neighboring points the tangent plane

		tr.clear();	
		for(unsigned int j = 0; j < neighbor_pts.size(); j ++){
			dVector vec = neighbor_pts[j] - pt;

			//cout<<vec<<endl;
			//cout<<neighbor_indices[j]<<endl;
			double x = CGAL::to_double(vec * tspace[0]);
			double y = CGAL::to_double(vec * tspace[1]);
			double z = CGAL::to_double(vec * tspace[2]);
			Vertex_handle_3d vh = tr.insert( Point_3d(x, y, z) );
			assert(vh != NULL);
			if(vh->get_origin_id() >= 0){
				//cout<<"warning: containing duplicated points"<<endl;
				fprintf(stderr,"Containing duplicated points\n");

			}
			vh->set_origin_id( neighbor_indices[j] );
		}

		map<unsigned int, double > vareas;
		for(Triangulation::Finite_vertices_iterator viter = tr.finite_vertices_begin(); viter != tr.finite_vertices_end(); viter ++){
			vareas[viter->get_origin_id()] = 0;
		}
		for(Triangulation::Finite_cells_iterator  citer = tr.finite_cells_begin(); citer != tr.finite_cells_end(); citer ++){
			Tetrahedron tetra( citer->vertex(0)->point(), citer->vertex(1)->point(), citer->vertex(2)->point(), citer->vertex(3)->point() );
			double carea = fabs( CGAL::to_double(tetra.volume()) );
			vareas[citer->vertex(0)->get_origin_id()] += carea / 4;
			vareas[citer->vertex(1)->get_origin_id()] += carea / 4;
			vareas[citer->vertex(2)->get_origin_id()] += carea / 4;
			vareas[citer->vertex(3)->get_origin_id()] += carea / 4;
		}
			
		double totalweight = 0;
		for(Triangulation::Finite_vertices_iterator viter = tr.finite_vertices_begin(); viter != tr.finite_vertices_end(); viter ++){
			unsigned int vid = viter->get_origin_id();
			if( vid == i ){
				assert( (viter->point() - Point_3d(0, 0, 0)) * (viter->point() - Point_3d(0, 0, 0)) < FT(1e-5) );
				continue;
			}

			double weight = exp(-CGAL::to_double( (viter->point() - Point_3d(0, 0, 0)) * (viter->point() - Point_3d(0, 0, 0))) / hh) / nmfactor;
			
			II.push_back(i + 1);
			JJ.push_back(vid + 1);
			SS.push_back(weight * vareas[vid]);
			//IIfout<<i + 1<<" ";
			//JJfout<<vid + 1<<" ";
			//SSfout<<weight * vareas[vid]<<" ";
			totalweight -= weight * vareas[vid];
		}
		
		II.push_back(i + 1);
		JJ.push_back(i + 1);
		SS.push_back(totalweight);
		//IIfout<<i + 1<<" ";
		//JJfout<<i + 1<<" ";
		//SSfout<<totalweight<<" ";
	}
	fprintf(stdout, "\n");
}

void generate_laplace_matrix_sparse_matlab_dimk(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS)
{
	//cout<<"generate_laplace_matrix_sparse_matlab_dimk"<<endl;
	fprintf(stdout, "Generate_laplace_matrix_sparse_matlab_dimk\n");

	typedef CGAL::Search_traits_d<KCd> Traits;
	typedef CGAL::Fair<Traits> Fair;
	typedef CGAL::Kd_tree<Traits, Fair> Tree;
	typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
	typedef KCd::FT FT;


	//const unsigned int tdim = 3;
	assert(tdim >= 1);

	//cout<<"h: "<<h<<" rho: "<<rho<<endl;
	double hh = h * h;
	double nmfactor = sqrt(M_PI * hh);
	for(unsigned int kk = 0; kk < tdim - 1; kk ++){
		nmfactor *= sqrt(M_PI * hh);
	}
	nmfactor *= (hh / 4);
	//cout<<"nmfactor: "<<nmfactor<<endl;
	fprintf(stdout, "nmfactor: %f\n", nmfactor);

	unsigned int np = pcloud.p_count();
	vector<dPoint> points;
	map<const dPoint, int, LessdPoint> pt2index;
	for(unsigned int i = 0; i < np; i ++){
		points.push_back( pcloud.point(i).coord() );
		pt2index.insert( make_pair(pcloud.point(i).coord(), i) );
	}
	Fair fair(10);
	Tree tree(points.begin(), points.end(), fair);
	points.clear();

	double eps = 1e-3 * h * rho;
	vector<dPoint> neighbor_pts;
	vector<int>	neighbor_indices;
	vector<dVector> tspace;
	Triangulation_d tr(tdim);

	for(unsigned int i = 0; i < np; i ++){

		//cout<<"i: "<<i<<"\r";
		fprintf(stdout, "i: %d\r", i);

		dPoint pt = pcloud.point(i).coord();
		
		//search the points within distance h * rho;
		neighbor_pts.clear();
		Fuzzy_sphere fs(pt, h * rho, FT(eps));
		tree.search(back_insert_iterator<vector<dPoint> >(neighbor_pts),  fs);

		//look up the indices of the neighboring points
		neighbor_indices.clear();

		bool find = false;
	   for(vector<dPoint>::iterator iter = neighbor_pts.begin(); iter != neighbor_pts.end(); iter ++){
  			//cout<<*iter<<endl;
			map<const dPoint, int, LessdPoint>::iterator iter_pt2index = pt2index.find(*iter);
			if( iter_pt2index != pt2index.end() ){
				neighbor_indices.push_back(iter_pt2index->second);
				if(iter_pt2index->second == (int)i){
					find = true;
				}
			}
			else{
				//cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
				fprintf(stdout, "Fail to find the index of the neighboring points.\n");
			}
   	}
		assert(neighbor_indices.size() >= 3 && find);

		//cout<<i<<"th neighbors: ";
		//for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
		//	cout<< neighbor_indices[j]<<" ";
		//}
		//cout<<endl;

		
		//estimate normal
		tspace.clear();
		est_tangent_space_lapack(pt, neighbor_pts, h, tdim, tspace);
		assert(tspace.size() == tdim);

		//cout<<"pt: "<<pt<<endl;
		//cout<<"tspace: "<<tspace[0]<<", "<<tspace[1]<<endl;

		tr.clear();
		map<dVertex_handle, ProjPoint> handle2vid;
		
		vector<double> coords;
		coords.resize(tdim, 0);
		//Point_d org(tdim, coords.begin(), coords.end());
		for(unsigned int j = 0; j < neighbor_pts.size(); j ++){
			dVector vec = neighbor_pts[j] - pt;

			//cout<<vec<<endl;
			//cout<<neighbor_indices[j]<<endl;
			for(unsigned int kk = 0; kk < tdim; kk ++){
				coords[kk] = CGAL::to_double(vec * tspace[kk]);
			}

			dVertex_handle vh = tr.insert( Point_d(tdim, coords.begin(), coords.end()) );
			assert(vh != NULL);
			//if(vh->get_origin_id() >= 0){
			//	cerr<<"warning: containing duplicated points"<<endl;
			//}
			//vh->set_origin_id( neighbor_indices[j] );
			map<dVertex_handle, ProjPoint>::iterator iter = handle2vid.find(vh);
			if(iter != handle2vid.end()){
				//cout<<"Error: Containing duplicated points"<<endl;
				fprintf(stdout, "Containing duplicated points\n");
			}
			else{
				//cout<<"dim: "<<tr.associated_point(vh).dimension()<<" "<<org.dimension()<<endl;
				//double d = CGAL::to_double( (vh->point() - org) * (vh->point() - org) );
				double d = CGAL::to_double( (tr.associated_point(vh) - CGAL::ORIGIN) * (tr.associated_point(vh) - CGAL::ORIGIN) );
				handle2vid.insert( make_pair(vh, ProjPoint(neighbor_indices[j], 0, d)) );
			}
		}

		for(Triangulation_d::Simplex_iterator siter = tr.simplices_begin(); siter != tr.simplices_end(); siter ++){
			vector<dPoint> points;
			for(unsigned int kk = 0; kk <= tdim; kk ++){
				points.push_back(tr.point_of_simplex(siter, kk));
			}
			double carea = fabs( simplex_sign_volume(points) );
			//double carea = 1;
			for(unsigned int kk = 0; kk <= tdim; kk ++){
				dVertex_handle vh = tr.vertex_of_simplex(siter, kk);
				map<dVertex_handle, ProjPoint>::iterator iter = handle2vid.find(vh);
				if(iter == handle2vid.end()){
					//cout<<"Error: can not find vertex"<<endl;
					fprintf(stdout, "Failed to find vertex\n");
				}
				else{
					iter->second = ProjPoint(iter->second.vid, iter->second.area + carea / (tdim + 1), iter->second.sqdist);
				}
			}
		}

		//--------------------------------------------------
		// //debug
		// cout<<"after computing area: "<<endl;
		// for(map<dVertex_handle, ProjPoint>::iterator iter = handle2vid.begin(); iter != handle2vid.end(); iter ++){
		// 	cout<< "coord: "<<iter->first->point() <<" vid: "<<iter->second.vid<<" area: "<<iter->second.area<<" dist: "<<iter->second.sqdist<<endl; 
		// }
		//-------------------------------------------------- 
		
		double totalweight = 0;
		
		for(map<dVertex_handle, ProjPoint>::iterator iter = handle2vid.begin(); iter != handle2vid.end(); iter ++){
			unsigned int vid = iter->second.vid;
			double area = iter->second.area;
			double sqdist = iter->second.sqdist;
			if( vid == i ){
				continue;
			}
			double weight = exp(-sqdist/hh) / nmfactor;
			
			//cout<<"vid: "<<vid<<" weight: "<<weight<<endl;
			II.push_back(i + 1);
			JJ.push_back(vid + 1);
			SS.push_back(weight * area);
			//IIfout<<i + 1<<" ";
			//JJfout<<vid + 1<<" ";
			//SSfout<<weight * area<<" ";
			totalweight -= weight * area;
		}

		II.push_back(i + 1);
		JJ.push_back(i + 1);
		SS.push_back(totalweight);
		//IIfout<<i + 1<<" ";
		//JJfout<<i + 1<<" ";
		//SSfout<<totalweight<<" ";
	}
	fprintf(stdout, "\n");
}


void generate_graph_laplace_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS)
{
	//cout<<"generate_graph_laplace_matrix_sparse_matlab"<<endl;
	fprintf(stdout, "generate_graph_laplace_matrix_sparse_matlab\n");

	typedef CGAL::Search_traits_d<KCd> Traits;
	typedef CGAL::Fair<Traits> Fair;
	typedef CGAL::Kd_tree<Traits, Fair> Tree;
	typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
	typedef KCd::FT FT;



	//cout<<"h: "<<h<<" rho: "<<rho<<endl;
	double hh = h * h;
	double nmfactor = sqrt(M_PI * hh);
	for(unsigned int kk = 0; kk < tdim - 1; kk ++){
		nmfactor *= sqrt(M_PI * hh);
	}
	nmfactor *= (hh / 4);
	//cout<<"nmfactor: "<<nmfactor<<endl;
	fprintf(stdout, "nmfactor: %f\n", nmfactor);

	unsigned int np = pcloud.p_count();
	vector<dPoint> points;
	map<const dPoint, int, LessdPoint> pt2index;
	for(unsigned int i = 0; i < np; i ++){
		points.push_back( pcloud.point(i).coord() );
		pt2index.insert( make_pair(pcloud.point(i).coord(), i) );
	}
	Fair fair(10);
	Tree tree(points.begin(), points.end(), fair);
	points.clear();

	double eps = 1e-3 * h * rho;
	vector<dPoint> neighbor_points;
	vector<int>	neighbor_indices;

	for(unsigned int i = 0; i < np; i ++){
		//cout<<"i: "<<i<<endl;
		fprintf(stdout, "i: %d\r", i);

		dPoint pt = pcloud.point(i).coord();
		//search the points within distance h * rho;
		neighbor_points.clear();
		Fuzzy_sphere fs(pt, h * rho, FT(eps));
		tree.search(back_insert_iterator<vector<dPoint> >(neighbor_points),  fs);

		//look up the indices of the neighboring points
		neighbor_indices.clear();

		bool find = false;
	   for(vector<dPoint>::iterator iter = neighbor_points.begin(); iter != neighbor_points.end(); iter ++){
  			//cout<<*iter<<endl;
			map<const dPoint, int, LessdPoint>::iterator iter_pt2index = pt2index.find(*iter);
			if( iter_pt2index != pt2index.end() ){
				neighbor_indices.push_back(iter_pt2index->second);
				if(iter_pt2index->second == (int)i){
					find = true;
				}
			}
			else{
				//cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
				fprintf(stdout, "Fail to find the index of the neighboring points\n");
			}
   	}
		assert(neighbor_indices.size() >= 3 && find);

		//cout<<i<<"th neighbors: ";
		//for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
		//	cout<< neighbor_indices[j]<<" ";
		//}
		//cout<<endl;

		double totalweight = 0;
		for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
			unsigned int pid = neighbor_indices[j];
			if(pid == i){
				continue;
			}

			dPoint pt_nb =  pcloud.point(pid).coord();
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) * ( 16.0 / (hh * hh * np));
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) * ( 4.0 / ( M_PI * hh * hh * np));
			double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) / (nmfactor * np);

			II.push_back(pid + 1);
			JJ.push_back(i + 1);
			SS.push_back(weight);
			//IIfout<<pid + 1<<" ";
			//JJfout<<i + 1<<" ";
			//SSfout<<weight<<" ";
		
			totalweight -= weight;
		}

		II.push_back(i + 1);
		JJ.push_back(i + 1);
		SS.push_back(totalweight);
		//IIfout<<i + 1<<" ";
		//JJfout<<i + 1<<" ";
		//SSfout<<totalweight<<" ";
	}
	fprintf(stdout, "\n");

	//--------------------------------------------------
	// IIfout<<"\n ";
	// JJfout<<"\n ";
	// SSfout<<"\n ";
	// IIfout.close();
	// JJfout.close();
	// SSfout.close();
	//-------------------------------------------------- 
}

void generate_arbdist_graph_laplace_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS)
{
	//cout<<"generate_graph_laplace_matrix_sparse_matlab"<<endl;
	fprintf(stdout, "generate_arbdist_graph_laplace_matrix_sparse_matlab\n");

	typedef CGAL::Search_traits_d<KCd> Traits;
	typedef CGAL::Fair<Traits> Fair;
	typedef CGAL::Kd_tree<Traits, Fair> Tree;
	typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
	typedef KCd::FT FT;



	//cout<<"h: "<<h<<" rho: "<<rho<<endl;
	double hh = h * h;
	double nmfactor = sqrt(M_PI * hh);
	for(unsigned int kk = 0; kk < tdim - 1; kk ++){
		nmfactor *= sqrt(M_PI * hh);
	}
	//nmfactor *= (hh / 4);
	//cout<<"nmfactor: "<<nmfactor<<endl;
	fprintf(stdout, "nmfactor: %f\n", nmfactor);

	unsigned int np = pcloud.p_count();
	vector<dPoint> points;
	map<const dPoint, int, LessdPoint> pt2index;
	for(unsigned int i = 0; i < np; i ++){
		points.push_back( pcloud.point(i).coord() );
		pt2index.insert( make_pair(pcloud.point(i).coord(), i) );
	}
	Fair fair(10);
	Tree tree(points.begin(), points.end(), fair);
	points.clear();

	double eps = 1e-3 * h * rho;
	vector<dPoint> neighbor_points;
	vector<int>	neighbor_indices;

	vector<double> density;
	for(unsigned int i = 0; i < np; i ++){
		density.push_back(0);
	}

	for(unsigned int i = 0; i < np; i ++){
		//cout<<"i: "<<i<<endl;
		fprintf(stdout, "i: %d\r", i);

		dPoint pt = pcloud.point(i).coord();
		//search the points within distance h * rho;
		neighbor_points.clear();
		Fuzzy_sphere fs(pt, h * rho, FT(eps));
		tree.search(back_insert_iterator<vector<dPoint> >(neighbor_points),  fs);

		//look up the indices of the neighboring points
		neighbor_indices.clear();

		bool find = false;
	   for(vector<dPoint>::iterator iter = neighbor_points.begin(); iter != neighbor_points.end(); iter ++){
  			//cout<<*iter<<endl;
			map<const dPoint, int, LessdPoint>::iterator iter_pt2index = pt2index.find(*iter);
			if( iter_pt2index != pt2index.end() ){
				neighbor_indices.push_back(iter_pt2index->second);
				if(iter_pt2index->second == (int)i){
					find = true;
				}
			}
			else{
				//cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
				fprintf(stdout, "Fail to find the index of the neighboring points\n");
			}
   	}
		assert(neighbor_indices.size() >= 3 && find);

		//cout<<i<<"th neighbors: ";
		//for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
		//	cout<< neighbor_indices[j]<<" ";
		//}
		//cout<<endl;

		//double totalweight = 0;
		for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
			unsigned int pid = neighbor_indices[j];
			if(pid == i){
				continue;
			}

			dPoint pt_nb =  pcloud.point(pid).coord();
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) * ( 16.0 / (hh * hh * np));
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) * ( 4.0 / ( M_PI * hh * hh * np));
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) / (nmfactor * np);
			double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) / (nmfactor);

			II.push_back(pid + 1);
			JJ.push_back(i + 1);
			SS.push_back(weight);
			//IIfout<<pid + 1<<" ";
			//JJfout<<i + 1<<" ";
			//SSfout<<weight<<" ";

			density[i] = density[i] + (weight / (np - 1));
			//totalweight -= weight;
		}
		
		//II.push_back(i + 1);
		//JJ.push_back(i + 1);
		//SS.push_back(totalweight);
		
		//IIfout<<i + 1<<" ";
		//JJfout<<i + 1<<" ";
		//SSfout<<totalweight<<" ";
	}
	fprintf(stdout, "\n");
	
	//--------------------------------------------------
	// for(unsigned int i = 0; i < np; i ++){
	// 	cout<<i<<" density "<<density[i]<<endl;
	// }
	//-------------------------------------------------- 

	//normalize
	assert(II.size() == JJ.size() && JJ.size() == SS.size());
	for(unsigned int i = 0; i < SS.size(); i ++){
		SS[i] = SS[i] / ( (np * (hh / 4)) * sqrt(density[II[i] - 1] * density[JJ[i] - 1]) );
	}
	//compute the diagonal elements
	

	vector<double> diags;
	for(unsigned int i = 0; i < np; i ++){
		diags.push_back(0);
	}
	for(unsigned int i = 0; i < SS.size(); i ++){
		//cout<<II[i]<<" ";
		diags[II[i] - 1] = diags[II[i] - 1] - SS[i];
	}

	for(unsigned int i = 0; i < np; i ++){
		II.push_back(i + 1);
		JJ.push_back(i + 1);
		SS.push_back(diags[i]);
	}

	//--------------------------------------------------
	// IIfout<<"\n ";
	// JJfout<<"\n ";
	// SSfout<<"\n ";
	// IIfout.close();
	// JJfout.close();
	// SSfout.close();
	//-------------------------------------------------- 
}


void generate_kernel_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS)
{
	fprintf(stdout, "generate_kernel_matrix_sparse_matlab\n");

	typedef CGAL::Search_traits_d<KCd> Traits;
	typedef CGAL::Fair<Traits> Fair;
	typedef CGAL::Kd_tree<Traits, Fair> Tree;
	typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
	typedef KCd::FT FT;



	//cout<<"h: "<<h<<" rho: "<<rho<<endl;
	double hh = h * h;

	//--------------------------------------------------
	// double nmfactor = sqrt(M_PI * hh);
	// for(unsigned int kk = 0; kk < tdim - 1; kk ++){
	// 	nmfactor *= sqrt(M_PI * hh);
	// }
	// nmfactor *= (hh / 4);
	// //cout<<"nmfactor: "<<nmfactor<<endl;
	// fprintf(stdout, "nmfactor: %f\n", nmfactor);
	//-------------------------------------------------- 

	unsigned int np = pcloud.p_count();
	vector<dPoint> points;
	map<const dPoint, int, LessdPoint> pt2index;
	for(unsigned int i = 0; i < np; i ++){
		points.push_back( pcloud.point(i).coord() );
		pt2index.insert( make_pair(pcloud.point(i).coord(), i) );
	}
	Fair fair(10);
	Tree tree(points.begin(), points.end(), fair);
	points.clear();

	double eps = 1e-3 * h * rho;
	vector<dPoint> neighbor_points;
	vector<int>	neighbor_indices;

	vector<double> density;
	density.resize(np, 0);

	for(unsigned int i = 0; i < np; i ++){
		//cout<<"i: "<<i<<endl;
		fprintf(stdout, "i: %d\r", i);

		dPoint pt = pcloud.point(i).coord();
		//search the points within distance h * rho;
		neighbor_points.clear();
		Fuzzy_sphere fs(pt, h * rho, FT(eps));
		tree.search(back_insert_iterator<vector<dPoint> >(neighbor_points),  fs);

		//look up the indices of the neighboring points
		neighbor_indices.clear();

		bool find = false;
	   for(vector<dPoint>::iterator iter = neighbor_points.begin(); iter != neighbor_points.end(); iter ++){
  			//cout<<*iter<<endl;
			map<const dPoint, int, LessdPoint>::iterator iter_pt2index = pt2index.find(*iter);
			if( iter_pt2index != pt2index.end() ){
				neighbor_indices.push_back(iter_pt2index->second);
				if(iter_pt2index->second == (int)i){
					find = true;
				}
			}
			else{
				//cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
				fprintf(stdout, "Fail to find the index of the neighboring points\n");
			}
   	}
		assert(neighbor_indices.size() >= 3 && find);

		//cout<<i<<"th neighbors: ";
		//for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
		//	cout<< neighbor_indices[j]<<" ";
		//}
		//cout<<endl;

		//double totalweight = 0;
		for(unsigned int j = 0; j < neighbor_indices.size(); j ++){
			unsigned int pid = neighbor_indices[j];

			dPoint pt_nb =  pcloud.point(pid).coord();
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) * ( 16.0 / (hh * hh * np));
			//double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh ) * ( 4.0 / ( M_PI * hh * hh * np));
			double weight = exp( -CGAL::to_double( (pt - pt_nb) * (pt - pt_nb) ) / hh );

			II.push_back(pid + 1);
			JJ.push_back(i + 1);
			SS.push_back(weight);
			
			density[i] += weight;

			//IIfout<<pid + 1<<" ";
			//JJfout<<i + 1<<" ";
			//SSfout<<weight<<" ";
			//totalweight -= weight;
		}
	}
	fprintf(stdout, "\n");

	assert(II.size() == JJ.size() && JJ.size() == SS.size());

	//normalize 1
	for(unsigned int i = 0; i < SS.size(); i ++){
		SS[i] /= (density[II[i] - 1] * density[JJ[i] - 1]) ;
	}
	
	//normalize 2
	density.clear();
	density.resize(np, 0);
	for(unsigned int i = 0; i < SS.size(); i ++){
		density[II[i] - 1] += SS[i];
	}
	for(unsigned int i = 0; i < SS.size(); i ++){
		SS[i] /= (density[II[i] - 1]);
	}

	
	//--------------------------------------------------
	// IIfout<<"\n ";
	// JJfout<<"\n ";
	// SSfout<<"\n ";
	// IIfout.close();
	// JJfout.close();
	// SSfout.close();
	//-------------------------------------------------- 
}


extern "C"  void dspev_(char &jobz, char &uplo, 
								long int &n, double *ap, 
								double *w, double *z, 
								long int &ldz, double *work, long int &info);
void est_tangent_space_lapack(const dPoint& pt, const vector<dPoint>& neighbor_pts, double h, unsigned int tdim, vector<dVector>& tspace)
{
	unsigned int dim = pt.dimension();
	assert(tdim <= dim);

	double* mat = new double[(dim + 1) * dim / 2];
	assert(mat != NULL);
	double* wr = new double[dim];
	assert(mat != NULL);
	double* vr = new double[dim * dim];
	assert(vr != NULL);

	double hh = h * h;
	memset(mat, 0, sizeof(double) * (dim + 1) * dim / 2);
	for(unsigned int i = 0; i < neighbor_pts.size(); i ++){
		dVector vec = neighbor_pts[i] - pt;
		double wt = exp( -CGAL::to_double(vec * vec) / hh );
		for(unsigned int j = 0; j < dim; j ++ ){
			for(unsigned int k = j; k < dim; k ++){
				mat[(k + 1) * k / 2 + j] += wt * vec[j] * vec[k];
			}
		}
	}

	// find out the eigenvalues and eigen vectors.
	char jobz = 'V', uplo = 'U';
	long int n = dim, ldz = dim, info = 0;
	double *work = new double[3 * dim];  
	assert(work != NULL);

	dspev_(jobz, uplo, n, mat, wr, vr, ldz, work, info);
	
	//cout<<"dspev_ info: "<<info<<endl;
	
	for(unsigned int i = 1; i <= tdim; i ++){
		tspace.push_back( dVector(dim, vr + (dim - i) * dim, vr + (dim - i + 1) * dim) );		
	}

	
	//cout<<"wr: ";
	//for(unsigned int j = 0; j < dim; j ++){
	//	cout<<wr[j]<<" ";
	//}
	//cout<<endl;
	
	delete []mat;
	delete []wr;
	delete []vr;
}

extern "C" void dgetrf_(long int &m, long int &n, 
							  double *a, long int &lda, 
							  long int* ipiv, long int &info);

double simplex_sign_volume(vector<dPoint> points)
{
	assert(points.size() > 0);
	unsigned int dim = points[0].dimension();

	assert(points.size() == dim + 1);

	double *A = new double[dim*dim];
	assert(A != NULL);
	for(unsigned int i = 0; i < dim; i ++){
		dVector vec = points[i+1] - points[0];
		for(unsigned int j = 0; j < dim; j ++){
			A[i * dim + j] = CGAL::to_double(vec[j]);
		}
	}
		
	long int m = dim, n = dim, lda = dim, info = 0;
	long int *ipiv = new long int[dim];
	assert(ipiv != NULL);

	dgetrf_(m, n, A, lda, ipiv, info);
	
	double volume = 1;
	double factor = 1;
	for(unsigned int i = 0; i < dim; i ++){
		volume *= A[i * dim + i];
		factor *= (i + 1);
	}
	volume /= factor;

	delete []ipiv;
	delete []A;
	
	//--------------------------------------------------
	// //debug
	// for(unsigned int i = 0; i <= dim; i ++){
	// 	cout<<points[i]<<endl;
	// }
	// cout<<"volume: "<<volume<<endl;
	// getchar();
	//-------------------------------------------------- 

	return volume;
}


