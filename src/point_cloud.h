#ifndef __POINT_CLOUD_H__
#define __POINT_CLOUD_H__

#include "datastructure.h"
//#include <mex.h>
using namespace std;

//-------------------------------------------------------------------
//VPCloud: vertex class for triangle mesh
//-------------------
#define VPCLOUD_FLAG_SELECTED	0x4


class VPCloud{

public:
  //Constructor
  	VPCloud(): _flag(0){}
  	VPCloud(const dPoint& co): _flag(0) {  _coord = co; }
  
  	VPCloud(const VPCloud& dP):_coord(dP.coord()), _flag(dP.flag()){}
//    	_coord = dP.coord(); 
//   	_flag = dP.flag();
//  	}
  	//Overloaded operator
  	VPCloud& operator = (const VPCloud& dP){
//		mexPrintf("This could be problem\n");
    	if (this != &dP){
			this->VPCloud::~VPCloud();
			new (this) VPCloud(dP);
		}
		return *this;
  	}
	//Function for accessing data
  	const dPoint& coord() const { return _coord; }
  	dPoint& coord() { return _coord; }
 	void set_coord(const dPoint& co) { _coord = co; }

  	bool check_flag(unsigned char f) { return _flag&f; }
  	void set_flag(unsigned char f) { _flag |= f; }
  	void un_set_flag(unsigned char f) { _flag &= (~f); }
  	unsigned char flag() const { return _flag; }
  	void copy_flag(unsigned char f) { _flag = f; }
  
private:
  	unsigned char _flag;
  	dPoint _coord;  
};

//-------------------------------------------------------------------


class PCloud{
public:
  	//Constructor
//  PCloud(){};
  	PCloud():_points(NULL), num_points(0){};
	PCloud(double *points, unsigned int np, unsigned int dim);
  	//Destructor
//  	virtual ~PCloud(){};
	~PCloud(){clear();}

  	//Function for accessing data
//  	unsigned int p_count() const { return _points.size(); }
    unsigned int p_count() const { return num_points; }

//  	unsigned int add_point(const VPCloud& v){ _points.push_back(v); return (_points.size() - 1);}
  	VPCloud& point(unsigned int i ) { return _points[i]; }
		
  	void copy_name(char* name) { sprintf(name, "%s", _name); }	
  	void set_name(char* name){ sprintf(_name, "%s", name );}
	
	unsigned int dd() const { return _d; }
	void set_dd(unsigned int d) { _d = d; }


	//Functions for output 
  	//void OutPCloudOffFile(char *filename);
  	//void OutPCloudOffFile(FILE *fp, double r, double g, double b, double a);
	void OutPCloud(char *filename);


  	//Functions for rendering
  	//void Render(float m[16]);
  	//void Render(float m[16], vector<double> fn, double min, double max, bool shownormal);
	//void Render_select_points(vector<unsigned int>& select_points);
	//void SelectRender_points();

  	//Functions for bounding box
  	void GetBBox();
  	dPoint pmin(){ return _pmin; }
  	dPoint pmax(){ return _pmax; }
  	double radius() {return sqrt( CGAL::to_double((_pmax - _pmin) * (_pmax - _pmin)) ) / 2; }

	//Average neighborhood size
  	double average_size(unsigned int k);
  	//Functions for Creation
 // 	bool ReadPointCloud(char *filename);
  	int ReadFromPCD(char *filename);

  	//Clear
//  	void clear(){ _points.clear();  }
	void clear(){ if (_points != NULL) delete [] _points;  }
            
private:
  unsigned int _d;
//  vector<VPCloud> _points;
  VPCloud *_points;
  int num_points;
  dPoint _pmin;
  dPoint _pmax;
  char _name[256];
};

#endif //__POINT_CLOUD_H__

