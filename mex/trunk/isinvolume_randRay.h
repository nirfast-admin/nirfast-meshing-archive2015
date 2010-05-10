#include "mex.h"
#include <vector>
#include <cmath>
#include <limits>
#include "vector.h"
#include <iostream>
#include <assert.h>
#include "geomath.h"

#ifndef ulong
#define ulong unsigned long
#endif

template<class T> inline T sqr(T x) { return (x)*(x); }
template<class T> inline T length(T x,T y,T z) { return sqrt(sqr(x)+sqr(y)+sqr(z)); }

struct points {
    double c[3];
};

unsigned int isinvolume_randRay(double *p0, const mxArray *mxP, const mxArray *mxT, double tiny, 
					   float (*facets_bbx)[6], int numPerturb, double minX, double maxX);
