// Usage: 
// [nodes] = tag_checkerbaord3d_mex(P, delta, llc, ds)

// To compile:
// On Windows platforms:
// mex -v -DWIN32 -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp
// mex -v -DWIN64 -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp

// On Apple Mac OSX:
// mex -v -DOSX -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp

// On Linux:
// mex -v -Dlinux -I./meshlib tag_checkerboard3d_mex.cpp meshlib/CStopWatch.cpp

#include "tag_checkerboard3d_mex.h"

#define P(i,j,k) P[(i)+nrow*(j)+ncol*nrow*(k)]
// #define row_state(i,j) row_state[(i)+npln*(j)]
#define row_state(i,j) row_state[npln*(i)+(j)]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	srand( (unsigned)time( NULL ) );
	CStopWatch mytimer;
	
	mwSize  nsubs; 
	const mwSize *dims;

	std::vector<mypoint> nodes;
	
	char *P = (char *) mxGetData(prhs[0]);
	double *delta = (double *) mxGetData(prhs[1]);
	double *llc = (double *) mxGetData(prhs[2]);
	double ds = *(mxGetPr(prhs[3]));
	
	#ifdef _mydebug
	mexPrintf(" Input data was read successfully!\n");
	#endif
	
	double dx, dy, dz;
	double xmin, ymin, zmin;
	int nrow, ncol, npln;
	
	dx = delta[0]; dy = delta[1]; dz = delta[2];
	xmin = llc[0]; ymin = llc[1]; zmin = llc[2];
	
	nsubs=mxGetNumberOfDimensions(prhs[0]);
	if (nsubs!=3)
		mexErrMsgTxt("Input checkerboard should be a 3D matrix!\n");
	
	dims = mxGetDimensions(prhs[0]);
	nrow = dims[0]; ncol = dims[1]; npln = (int) dims[2];
	ULONG no_tagged_pixels=0;
	ULONG cc=0;
	
	#ifdef _mydebug
	mexPrintf(" Read dimensions successfully!\n");
	#endif

	mytimer.startTimer();
	
	#ifdef _mydebug
	mexPrintf(" Timer initialized successfully!\n");
	#endif
	
	#ifdef _mydebug
	mexPrintf(" nrow = %ld , npln = %ld\n",nrow,npln);
	#endif
	int *row_state = new int[nrow*npln];
	for (int i=0; i<nrow; ++i) {
		for (int j=0; j<npln; ++j) {
	#ifdef _mydebug
			mexPrintf(" i = %ld , j = %ld\n", i, j);
	#endif
			row_state(i,j) = 0;
		}
	}
	#ifdef _mydebug
	mexPrintf(" Initialized row_state() successfully!\n");
	#endif

	mexPrintf("  Tagging Possible Interior Nodes...");
	
	while (true) {
		int ii, jj, kk;
        ii = myrand(nrow);
        kk = myrand(npln);
		jj = row_state(ii,kk);
        
		if (jj >= ncol) {
			++cc;
			continue;
		}
		
		int idx = jj;
		for (idx=jj; idx<ncol; ++idx) {
			if (P(ii,idx,kk) == 0) {
				P(ii,idx,kk) = node_code;
				
				mypoint foo;
				foo.coords[0] = (idx-2)*dx+xmin;
				foo.coords[1] = (nrow-ii-1)*dy+ymin;
				foo.coords[2] = (npln-kk-1)*dz+zmin;
				nodes.push_back(foo);
				++no_tagged_pixels;
                
				double dd = ds;
				int  dI = round(dd/dx) - 1;
				for (int i=_stupidMS_max(0,ii-dI); i<_stupidMS_min(nrow,ii+dI); ++i) {
					for (int j=_stupidMS_max(0,idx-dI); j<_stupidMS_min(ncol, idx+dI); ++j) {
						for (int k=_stupidMS_max(0,kk-dI); k<_stupidMS_min(npln,kk+dI); ++k) {
							if (P(i,j,k)==0)
								P(i,j,k) = NA;
						}
					}
				}
				break;
			}
		}
		row_state(ii,kk) = idx + 1;
		if (ReadyForTermination(row_state,nrow,ncol,npln))
			break;
	}
	mytimer.stopTimer();
	mexPrintf(" done!\n");
	
	#ifdef _mydebug
	mexPrintf(" Tagged the pixels successfully\n");
	#endif

	if (no_tagged_pixels==0)
		mexErrMsgTxt("checkerboard3d:EmptyMatrix','Based on input matrix, checkerboard3d can not deploy any nodes!");
	
	/*mexPrintf("  Number of calling rand() function redundantly: %d\n",cc);
	mexPrintf("  nrow*ncol*npln = %d\n", nrow*ncol*npln);
	mexPrintf("  Time elapsed in tagging : %.6g\n", mytimer.getElapsedTime());*/
	
	plhs[0] = mxCreateNumericMatrix(nodes.size(), 3, mxDOUBLE_CLASS, mxREAL);
	double *foo = (double *) mxGetData(plhs[0]);
	for (ULONG i=0; i<nodes.size(); ++i) {
		for (int j=0; j<3; ++j) {
			foo[i+nodes.size()*j] = nodes[i].coords[j];
		}
	}
	
	#ifdef _mydebug
	mexPrintf(" Assigned ouput successfully\n");
	#endif

	nodes.clear();
    delete [] row_state;
	
}


	
bool ReadyForTermination(int *row_state, int& nrow, int& ncol, int& npln) {
	for (int i=0; i<nrow; ++i) {
		for (int j=0; j<npln; ++j) {
			if (row_state(i,j) < ncol)
				return false;
		}
	}
	return true;
}
