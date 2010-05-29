// Matlab usage:
// st = PointInPolyhedron_mex(P, polygons, points, eps)
// st = 0 : Outside
// st = 1 : Inside
// st = 2 : On the boundary of polyhedron
// Input:
// P = [x y z]
// polygons : node connectivy 
//            n1 n2 n3 ...
//            n4 n5 n2 ...
// points   : coordiantes of polyhedron vertices.

// To compile this in Matlab use:

// Windows:
// mex -v -I./header -DCPU86 -DWIN32 PointInPolyhedron_mex.cpp header/polyhedron2BSP.cpp header/CPoint.cpp header/CVector.cpp header/Plane3D.cpp header/BSPNode.cpp header/MeshIO.cpp header/FileOperation.cpp header/CPolygon.cpp header/predicates.cpp header/CStopWatch.cpp

// Linux:
// mex -v -I./header -DLINUX PointInPolyhedron_mex.cpp header/polyhedron2BSP.cpp header/CPoint.cpp header/CVector.cpp header/Plane3D.cpp header/BSPNode.cpp header/MeshIO.cpp header/FileOperation.cpp header/CPolygon.cpp header/predicates.cpp header/CStopWatch.cpp

// Mac OSX
// mex -v -I./header         PointInPolyhedron_mex.cpp header/polyhedron2BSP.cpp header/CPoint.cpp header/CVector.cpp header/Plane3D.cpp header/BSPNode.cpp header/MeshIO.cpp header/FileOperation.cpp header/CPolygon.cpp header/predicates.cpp header/CStopWatch.cpp


#include "PointInPolyhedron_mex.h"
#include "init.h"

#define points(i,j) points[(i)+np*(j)]
#define ele(i,j) ele[(i)+ne*(j)]
#define t(i,j) t[(i)+ne*(j)]


int PointInSolidSpace(BSPNode *node, Point& p, double PlaneTHK = TinyZero);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    CheckArgsIn(nlhs, plhs, nrhs, prhs);

	
    CStopWatch mytimer;
    bool debug = false;
    
    unsigned long ne = mxGetM(prhs[1]);
    int nnpe = mxGetN(prhs[1]);
	double *t;
    ULONG *ele = new ULONG[ne*nnpe];
    if (mxIsDouble(prhs[1])) {
        t = (double *) mxGetData(prhs[1]);
        for (ULONG i=0; i<ne; ++i) {
            for (int j=0; j<nnpe; ++j)
                ele(i,j) = (ULONG) t(i,j);
        }
    }
    else
        mexErrMsgTxt("PointInPolyhedron_mex: input ele list should be in double\n");

	double *points = (double *) mxGetData(prhs[2]);

	
	unsigned long np = mxGetM(prhs[2]);

	double eps = 0.;
	if (nrhs>=4) {
		eps = *(mxGetPr(prhs[3]));
	}
	else
		eps = TinyZero;
    
    //mexPrintf("TinyZero is : %.12f\n",TinyZero);

	double *p = (double *) mxGetData(prhs[0]);
	ULONG nqueryp = mxGetM(prhs[0]);

    if (debug)
        mexPrintf("Creating the output array\n");
	plhs[0] = mxCreateNumericMatrix(nqueryp, 1, mxINT8_CLASS, mxREAL);
	char *c = (char *) mxGetData(plhs[0]);

	Polyhedron2BSP myconverter;
    if (debug)
        mexPrintf("Initializing the polyhedron converter\n");
	myconverter.InitFromMatlabMex(points, ele, np, ne, nnpe);
	
    mexPrintf("\nBuilding the BSP tree...");
    
	mytimer.startTimer();
	double foothk = myconverter.GetPlaneThickness();
    //mexPrintf("\nfoothk: %10.16g\n",foothk);
    if (debug)
        mexPrintf("Runnig the query...\n");
	for (ULONG i=0; i<nqueryp; ++i) {
		Point qfoo(p[i], p[i + nqueryp], p[i + nqueryp + nqueryp]);
        *(c+i) = myconverter.IsInside(qfoo, foothk);
		// *(c+i) = PointInSolidSpace(root, qfoo, foothk);
	}
	mytimer.stopTimer();
    mexPrintf(" done! (Time: %.12f sec.)\n",mytimer.getElapsedTime());

	// Testing to see if all polygons of the input polyhedron have been used as a leaf in the BSP tree
	ULONG not_visited = 0;
	for (ULONG i=0; i<myconverter.polygonmarker.size(); ++i) {
		if (myconverter.polygonmarker[i] == false) {
			++not_visited;
		}
	}
	//mexPrintf("\n(%d out of %d polygons were not used in the BSP tree construction!)\n",not_visited,(ULONG)myconverter.polygonmarker.size());


	//mexPrintf("\n  Releasing tree memory");
	myconverter.DeleteTree();
	//mexPrintf("... done.\n\n");
}




bool CheckArgsIn(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if (nrhs < 3 || nrhs > 4) {
        mexPrintf("nargin = %d\n",nrhs);
        mexErrMsgTxt("PointInPolyhedron_mex: You need to enter 4 input args!");
    }
    if (nlhs!=1) {
        mexErrMsgTxt("PointInPolyhedron_mex: Output argument is not specified!");
    }
	if (mxGetN(prhs[2])!=3) {
		mexErrMsgTxt("PointInPolyhedron_mex can only be used for 3D shells!\n");
	}
    bool f1=true;
    for (int i=0; i<nrhs; f1=f1 && !mxIsComplex(prhs[i]), ++i);
    if (!f1) {
        mexErrMsgTxt("PointInPolyhedron_mex: Input arguments are complex!");
    }
    for (int i=0; i<nrhs; f1=f1 && (!mxIsEmpty(prhs[i]) || i==5), ++i);
    if (!f1) {
        mexErrMsgTxt("PointInPolyhedron_mex: Input arguments are empty!");
    }
    f1=true;
    for (int i=0; i<nrhs; f1=f1 && !mxIsChar(prhs[i]), ++i);
    if (!f1) {
        mexErrMsgTxt("PointInPolyhedron_mex: Input arguments are string!");
    }
	
	if (/* P */!mxIsDouble(prhs[0]) || (nrhs==4 && /*eps*/!mxIsDouble(prhs[3])) )
        mexErrMsgTxt("P and eps need to be double");

	if (!mxGetPr(prhs[2])) {
		mexErrMsgTxt("points need to be double");
	}
 
	return true;
}

