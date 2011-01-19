#include "involume_mex.h"

#define p(i,j) p[(i)+np*(j)]
#define ele(i,j) ele[(i)+ne*(j)]
#define qp(i,j) qp[(i)+nqp*(j)]

/* To compile this file use:
 * Compile without OpenMP support:
 *-------------
* For Windows:
* mex -v -DWIN32 -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp
For Linux/Mac:
 * mex -v -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp
 *
 *Compile with OpenMP support:
 *For Windows:
 *mex -v COMPFLAGS="$COMPFLAGS /openmp" LINKFLAGS="$LINKFLAGS /openmp" -DWIN32 -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp
 *
 *For Linux/Mac:
 *mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I./meshlib involume_mex.cpp isinvolume_randRay.cpp meshlib/geomath.cpp meshlib/vector.cpp
 */

#ifdef _OPENMP
#include "omp.h"
#endif

// st = involume_mex(qp, ele, p, ntries, facets_bbx, xMin, xMax, tiny)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs!=8 || nlhs!=1) {
        mexPrintf("st = involume_mex(qp, ele, p, facets_bbx, ntries, xMin, xMax, tiny)");
        mexPrintf("It needs 8 inputs and 1 output\n");
        mexErrMsgTxt("involume_mex: Terminatin...");
    }
    if (!mxIsDouble(_innode))
        mexErrMsgTxt("involume_mex: input argument 2 should be of 'double' type");
    
    bool debug=false;
    
    if (mxGetN(_inele)!=3 || mxGetN(_inqp) != 3)
        mexErrMsgTxt("involume_mex: is intended for 3D!");
    double *ele;
    if (mxIsDouble(_inele))
        ele = (double *) mxGetData(prhs[1]);
    else
        mexErrMsgTxt("orient_surface: t needs to be 'double'!");
    
	/*ULONG ne = mxGetM(_inele);
    ULONG np = mxGetM(_innode);*/
    ULONG nqp = (ULONG) mxGetM(_inqp);
    //double *p = mxGetPr(_innode);
    double *qp = mxGetPr(_inqp);

    double tmpqp[3] = {0., 0., 0.};
    
	double ntries = mxGetScalar(_inntries);
    double xMin = mxGetScalar(_inxmin);
    double xMax = mxGetScalar(_inxmax);
    double tiny = mxGetScalar(_intiny);
    
    //     check to see if we need to create facets_bbx
	double *facets_bbx;
   
    if (debug)
        mexPrintf("Getting facets_bbx\n");
	if (!mxIsEmpty(_infacets))
		facets_bbx = (double *) mxGetData(_infacets);
	else
		mexErrMsgTxt("involume_mex: facets_bbx is empty!\n");
	
    _outst = mxCreateNumericMatrix(nqp, 1, mxUINT8_CLASS, mxREAL);
    unsigned char *st = (unsigned char *) mxGetData(_outst);
    
    if (debug)
        mexPrintf("Entering the loop to call isinvolume_ranRay\n");
	
#ifdef _OPENMP	
	omp_set_num_threads(omp_get_num_procs());
    if (debug)
        mexPrintf("\n    CPUs Available: %d\n\n",omp_get_num_procs());
#endif
	long i;
	int j;
    std::vector<ULONG> int_facets;
    std::vector<points> int_points;
#ifdef _OPENMP
	#pragma omp parallel default(none) \
	    private(i,j) \
        firstprivate(int_facets,int_points,tmpqp) \
	    shared(st,qp,nqp,prhs,tiny,facets_bbx,ntries,xMin,xMax)
#endif
	{
#ifdef _OPENMP
        #pragma omp for nowait
#endif
        for (i=0; i<nqp; ++i) {
    	    for (j=0; j<3; tmpqp[j] = qp(i,j), ++j);
    		st[i] = isinvolume_randRay(tmpqp,prhs[2],prhs[1],tiny,facets_bbx,(int)ntries,xMin,xMax,int_facets,int_points);
    	}
	}
}







