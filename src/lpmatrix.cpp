//#include <mex.h>
#include "comp_llpmatrix.h"
#include "lpmatrix.h"


void generate_pcdlaplace_matrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	PCloud pcloud(points, np, dim);

	//--------------------------------------------------
   //ofstream fout; 
	//fout.open("points");
	//for(unsigned int i = 0; i < np; i ++){ 
   //	for(unsigned int j = 0; j < dim; j ++){ 
	//  	fout<< points[j * np + i] <<" ";
	//	}
	//	fout<<endl;
	//}
	//pcloud.OutPCloud("pcd");
	//--------------------------------------------------
	double h;
	double avers = pcloud.average_size(nn);
	printf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	printf("h: %f\n", h);

	if(tdim == 2){
		generate_laplace_matrix_sparse_matlab_dim2(pcloud, h, rho, IIV, JJV, SSV);
		//generate_laplace_matrix_sparse_matlab_dimk(pcloud, h, rho, tdim, IIV, JJV, SSV);
	}
	else if(tdim == 3){
		generate_laplace_matrix_sparse_matlab_dim3(pcloud, h, rho, IIV, JJV, SSV);
		//generate_laplace_matrix_sparse_matlab_dimk(pcloud, h, rho, tdim, IIV, JJV, SSV);
	}
	else{
		generate_laplace_matrix_sparse_matlab_dimk(pcloud, h, rho, tdim, IIV, JJV, SSV);
	}

}

void generate_graphlaplace_matrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	PCloud pcloud(points, np, dim);

	//--------------------------------------------------
   //ofstream fout; 
	//fout.open("points");
	//for(unsigned int i = 0; i < np; i ++){ 
   //	for(unsigned int j = 0; j < dim; j ++){ 
	//  	fout<< points[j * np + i] <<" ";
	//	}
	//	fout<<endl;
	//}
	//pcloud.OutPCloud("pcd");
	//--------------------------------------------------
	double h;
	double avers = pcloud.average_size(nn);
	printf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	printf("h: %f\n", h);
	generate_graph_laplace_matrix_sparse_matlab(pcloud, h, rho, tdim, IIV, JJV, SSV);
}

void generate_kernelmatrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV, double& h)
{
	PCloud pcloud(points, np, dim);

	//--------------------------------------------------
   //ofstream fout; 
	//fout.open("points");
	//for(unsigned int i = 0; i < np; i ++){ 
   //	for(unsigned int j = 0; j < dim; j ++){ 
	//  	fout<< points[j * np + i] <<" ";
	//	}
	//	fout<<endl;
	//}
	//pcloud.OutPCloud("pcd");
	//--------------------------------------------------
	double avers = pcloud.average_size(nn);
	printf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	printf("h: %f\n", h);
	generate_kernel_matrix_sparse_matlab(pcloud, h, rho, tdim, IIV, JJV, SSV);
}

void generate_arbdistgraphlaplace_matrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	PCloud pcloud(points, np, dim);

	//--------------------------------------------------
   //ofstream fout; 
	//fout.open("points");
	//for(unsigned int i = 0; i < np; i ++){ 
   //	for(unsigned int j = 0; j < dim; j ++){ 
	//  	fout<< points[j * np + i] <<" ";
	//	}
	//	fout<<endl;
	//}
	//pcloud.OutPCloud("pcd");
	//--------------------------------------------------
	double h;
	double avers = pcloud.average_size(nn);
	printf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	printf("h: %f\n", h);
	generate_arbdist_graph_laplace_matrix_sparse_matlab(pcloud, h, rho, tdim, IIV, JJV, SSV);
}


