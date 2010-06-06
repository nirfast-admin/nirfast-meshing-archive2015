/*
 *  intersect_ray_shell_mex.h
 *  
 *
 *  Created by Hamid Ghadyani on 12/18/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __intersect_ray_shell_mex_h
#define __intersect_ray_shell_mex_h

#ifndef sign(x)
#define sign(x) (fabs(x)/(x))
#endif
#define ulong unsigned long int

#include "mex.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <assert.h>
#include "vector.h"

//#include "segseg.h"
//#include "ray_triangle_coplanar.h"
#include "geomath.h"
//#include "intersect_ray_triangle_moller.h"

#ifndef TinyZero
#define TinyZero 1e-10
#endif 
struct points {
    double c[3];
};
unsigned long m_indexing;


bool CheckArgsIn(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//unsigned long int *GetEdgeTriConnList(mxArray *list);
ulong GetNeighborTriangle(ulong edge[2], ulong tri, const mxArray *list, ulong *indexing);
bool IsValidCrossing(ulong idx, int first_st, std::vector<std::vector<double> > first_intpnts, 
					 double *rp1, double *rp2, double *mydir, double tiny,
					 std::vector<points> &int_pnts, std::vector<int> &int_status, std::vector<bool> &int_facets,
					 double *shell_normals, const mxArray *list, ulong *indexing,
					 double *t, double *p, ulong ne, ulong np);

#endif

