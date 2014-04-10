#ifndef __LPMATRIX_H__
#define __LPMATRIX_H__

#include <vector>
using namespace std;

void generate_pcdlaplace_matrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV);

void generate_graphlaplace_matrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV);

void generate_kernelmatrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV, double& h);

void generate_arbdistgraphlaplace_matrix_sparse_matlab(double *points, unsigned int np, unsigned int dim, unsigned int tdim, unsigned int htype, unsigned int nn, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV);

#endif //__LPMATRIX_H__
