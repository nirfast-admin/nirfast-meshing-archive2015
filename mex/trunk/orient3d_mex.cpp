#include "mex.h"
#include "geomath.h"
#ifndef TinyZero
#define TinyZero 1e-10
#endif 

// inflag = orient3d_mex(pa, pb, pc, pd)
/* Return a positive value if the point pd lies below the      */
/* plane passing through pa, pb, and pc; "below" is defined so */
/* that pa, pb, and pc appear in counterclockwise order when   */
/* viewed from above the plane.  Returns a negative value if   */
/* pd lies above the plane.  Returns zero if the points are    */
/* coplanar.  The result is also a rough approximation of six  */
/* times the signed volume of the tetrahedron defined by the   */
/* four points.                                                */

// To mex under all platforms:
// mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" orient3d_mex.cpp meshlib/geomath.cpp meshlib/vector.cpp meshlib/predicates.cpp -I./meshlib

int orient3d(double *pa, double *pb, double *pc, double *pd);

#ifdef _OPENMP
#include "omp.h"
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs!=1 || nrhs!=4)
        mexErrMsgTxt("orient3d_mex: needs 4 input and 1 output\n");
    
    unsigned long nvert = mxGetM(prhs[0]);
    if (nvert < 1)
      return;
    if (mxGetM(prhs[0]) != mxGetM(prhs[1]) ||
        mxGetM(prhs[0]) != mxGetM(prhs[2]) ||
        mxGetM(prhs[0]) != mxGetM(prhs[3]) ||
        mxGetM(prhs[1]) != mxGetM(prhs[2]) ||
        mxGetM(prhs[1]) != mxGetM(prhs[3]) ||
        mxGetM(prhs[2]) != mxGetM(prhs[3]))
      mexErrMsgTxt("orient3d_mex: input points should be the same size!\n");
    if (mxGetN(prhs[0]) != 3 || mxGetN(prhs[1]) != 3 ||
        mxGetN(prhs[2]) != 3 || mxGetN(prhs[3]) != 3 )
      mexErrMsgTxt("orient3d_mex: need 3D input vertices!\n");

    plhs[0] = mxCreateNumericMatrix(nvert, 1, mxINT8_CLASS, mxREAL);
    char *c = (char *) mxGetData(plhs[0]);
    long i = 0, j = 0;
    double foo = 0.;
    double pa[3], pb[3], pc[3], pd[3];
    for (j=0; j<3; ++j) {
      pa[j] = 0.;
      pb[j] = 0.;
      pc[j] = 0.;
      pd[j] = 0.;
    }

    #ifdef _OPENMP
    int nt = omp_get_num_procs();
    #endif
    #pragma omp parallel default(none) \
    firstprivate(i, j, foo, pa, pb, pc, pd) \
    shared (nvert, prhs, plhs, c, nt)
    {
      #pragma omp master
      mexPrintf("Number of threads: %d\n", nt);

      #pragma omp for
      for (i=0; i<nvert; ++i) {
        for (j=0; j<3; ++j) {
          pa[j] = mxGetPr(prhs[0])[i + j*nvert];
          pb[j] = mxGetPr(prhs[1])[i + j*nvert];
          pc[j] = mxGetPr(prhs[2])[i + j*nvert];
          pd[j] = mxGetPr(prhs[3])[i + j*nvert];
        }
        foo = orient3d(pa, pb, pc, pd);
        if (IsEqual(foo, 0., TinyZero))
          *(c+i) = 0;
        else if (foo > 0)
          *(c+i) = 1;
        else
          *(c+i) = -1;
      }
    }
}
