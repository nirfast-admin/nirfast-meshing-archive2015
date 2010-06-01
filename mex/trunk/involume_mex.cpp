#include "involume_mex.h"

#define p(i,j) p[(i)+np*(j)]
#define ele(i,j) ele[(i)+ne*(j)]
#define qp(i,j) qp[(i)+nqp*(j)]
/* To compile this file use:
* For Windows:
* mex -v -DWIN32 -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp
For Linux/Mac:
 * mex -v -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp
 */

// st = involume_mex(qp, ele, p, ntries, xMin, xMax, tiny)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs!=7 || nlhs!=1) {
        mexErrMsgTxt("st = involume_mex(qp, ele, p, facets_bbx, ntries, xMin, xMax, tiny)");
        mexErrMsgTxt("involume_mex: Terminatin...");
    }
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("involume_mex: input argument 2 should be of 'double' type");
    
    bool debug=false;
    
    int dim = mxGetN(prhs[1]);
    if (mxGetN(prhs[1])!=3 || mxGetN(prhs[0]) != 3)
        mexErrMsgTxt("involume_mex: is intended for 3D!");
    double *ele;
    if (mxIsDouble(prhs[1]))
        ele = (double *) mxGetData(prhs[1]);
    else
        mexErrMsgTxt("orient_surface: t needs to be 'double'!");
    
	ulong ne = mxGetM(prhs[1]);
    ulong np = mxGetM(prhs[2]);
    ulong nqp = mxGetM(prhs[0]);
    double *p = mxGetPr(prhs[2]);
    double *qp = mxGetPr(prhs[0]);
    double tmpqp[3] = {0.,0.,0.};
    
	double ntries = mxGetScalar(prhs[3]);
    double xMin = mxGetScalar(prhs[4]);
    double xMax = mxGetScalar(prhs[5]);
    double tiny = mxGetScalar(prhs[6]);
    
    //     check to see if we need to create facets_bbx
    float (*facets_bbx)[6];
   
    if (debug)
        mexPrintf("Calculating facets_bbx\n");
    facets_bbx = new float[ne][6];
    for (ulong i=0; i<ne; ++i) {
        ulong n1 = (ulong) ele(i,0);
        ulong n2 = (ulong) ele(i,1);
        ulong n3 = (ulong) ele(i,2);
        double tmp;
        for (ulong j=0; j<3; ++j) {
            tmp = std::min(p(n1-1,j),p(n2-1,j));
            facets_bbx[i][j]=std::min(tmp,p(n3-1,j));
            tmp = std::max(p(n1-1,j),p(n2-1,j));
            facets_bbx[i][j+3]=std::max(tmp,p(n3-1,j));
        }
    }
    
    plhs[0] = mxCreateNumericMatrix(nqp, 1, mxUINT8_CLASS, mxREAL);
    unsigned char *st = (unsigned char *) mxGetData(plhs[0]);
    
    if (debug)
        mexPrintf("Entering the loop to call isinvolume_ranRay\n");
    for (ulong i=0; i<nqp; ++i) {
        for (int j=0; j<3; tmpqp[j] = qp(i,j), ++j);
        st[i] = isinvolume_randRay(tmpqp,prhs[2],prhs[1],tiny,facets_bbx,ntries,xMin,xMax);
    }
}

