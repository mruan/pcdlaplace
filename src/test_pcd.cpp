//#include <mex.h>
#include <math.h>
#include <vector>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>

#include "point_cloud.h"

int main(int argc, char* argv[])
{
  
  PCloud pcloud;
  printf("Try to load file %s\n", argv[1]);
  int np = pcloud.ReadFromPCD(argv[1]);//np < 0 if error 

  int tdim = 2, dim = 3;	
  unsigned int nn = 10, htype = 0;
  double hs = 2, rho = 3;

  printf("np: %d, dim: %d, tdim: %d, ", np, dim, tdim);
  printf("htype: %d, nn: %d, hs: %.2f, rho: %.2f, ", htype, nn, hs, rho);
  return 0;
}
