#include "isinvolume_randRay.h"

#ifndef p
#define p(i,j) p[(i)+np*(j)]
#define t(i,j) t[(i)+ne*(j)]
#endif

#ifndef _SIGN
#define _SIGN(x) (fabs(x)/(x))
#endif


// This routine should be called from a Matlab mex function ONLY!
unsigned char isinvolume_randRay(double *p0, const mxArray *mxP, const mxArray *mxT, double tiny, 
					   float (*facets_bbx)[6], int numPerturb, double minX, double maxX)
{
   
    unsigned char st;
    bool debug=false;
        
    // get p0 and list of shell points
    if (/*rp1 & rp2 */!mxIsDouble(mxP))
        mexErrMsgTxt("isinvolume_randRay: p0, p and numPerturb need to be of type 'double'!");
	double *p = mxGetPr(mxP);
    if (debug)
        mexPrintf("p0 = %lf %lf %lf\n",p0[0],p0[1],p0[2]);
    // get list of shell elements
    double *t;
    if (mxIsDouble(mxT))
        t = (double *) mxGetData(mxT);
    else
        mexErrMsgTxt("isinvolume_randRay: t needs to be 'double'!");
    
    // extract number of points and elements
    unsigned long np=mxGetM(mxP);
	unsigned long ne=mxGetM(mxT);


    // c++ vector to hold growing list of intersection points
	std::vector<points> intpnts;
    points p1;
    
    if (debug)
        mexPrintf("Calculating ray's bbx\n");
    
    // mxArrays to be passed down to intersect_RayTriangle
    // double *pr;
    if (debug)
        mexPrintf("Duplicating input arguments\n");
    
    
    double rp1[3], rp2[3];
    for (int i=0; i<3; rp1[i]=p0[i], rp2[i]=p0[i], ++i);
    // and add 2*AbsMaxX (as infinite point)
	rp2[0] = fabs(maxX-minX)*10 + maxX;
    
    unsigned long int nintpnts;
	bool forever = true;
    bool brokenloop = false;
    bool rayPerturbed=false;
    double R;
    
    if (debug)
        mexPrintf("Entering the main loop\n");
	while (forever) {
        for (ulong i=0; i<ne; ++i) {
            if (!rayPerturbed) {
                double MaxX = facets_bbx[i][3]; 
                double MaxY = facets_bbx[i][4]; double MinY = facets_bbx[i][1];
                double MaxZ = facets_bbx[i][5]; double MinZ = facets_bbx[i][2];
                if (debug) {
//                     mexPrintf("Ray Is Perturbed! i = %d\n",i);
                }
                if (p0[0]>MaxX || p0[1]>MaxY || p0[1]<MinY || p0[2]>MaxZ || p0[2]<MinZ)
                    continue;
            }
            // get coordinates of tp1, tp2 and tp3 into mxArrays of *rhs1[]
            double tp1[3], tp2[3], tp3[3];
            for (int k=0; k<3; ++k) {
		ulong idx = (ulong) t(i,0) - 1;
		assert(idx>=0 && idx<=np);
                tp1[k] = p(idx, k);
		idx = (ulong) t(i,1) - 1;
		assert(idx>=0 && idx<=np);
                tp2[k] = p(idx, k);
		idx = (ulong) t(i,2) - 1;
		assert(idx>=0 && idx<=np);
                tp3[k] = p(idx, k);
            }

            if (debug) {
//                 mexPrintf("Calling intersect_RayTriangle\n");
            }
            // [st I]=intersect_RayTriangle(rp1, rp2, tp1, tp2, tp3, global_tiny)
            double ipnt[3];
            int ret = 
            intersect_RayTriangle(rp1, rp2, tp1, tp2, tp3, ipnt, tiny);
            if (ret==1) {
                for (int j=0;j<3;p1.c[j]=ipnt[j],++j);
				intpnts.push_back(p1);
                if (debug)
                    mexPrintf("i=%d \t intp = %lf %lf %lf\n",i,p1.c[0], p1.c[1], p1.c[2]);
                continue;
            }
            else if (ret==10 || ret==11 || ret==12 || ret==20 || ret==21 || ret==22
                     || ret>100) {
                if (!rayPerturbed) {
                    R = -std::numeric_limits<double>::max();
                    for (ulong j=0;j<np;++j) {
                        double temp = length(p0[0]-p(j,0), p0[1]-p(j,1), p0[2]-p(j,2));
                        if (temp > R)
                            R = temp;
                    }
                    R *= 100;
                    rayPerturbed=true;
                    if (debug)
                        mexPrintf("Calculated R is : %lf\n",R);
                }
                brokenloop=true;
				intpnts.clear();
                break;
            }
        }
        nintpnts = intpnts.size();
        if (!brokenloop) {
            (nintpnts%2==1) ? st=1 : st=0;
            if (debug) {
                mexPrintf("numPerturb is %d\n",numPerturb);
                mexPrintf("nintpnts=%d st=%d\n",nintpnts,st);
            }
            forever=false;
        }
        else if (brokenloop && numPerturb>0) {
            double u[3] = { rand(), rand(), rand() };
            double temp = length(u[0],u[1],u[2]);
            for (int j=0; j<3; u[j] /= temp/*, mexPrintf("u= %.12f\n",u[j])*/, ++j);
            for (int j=0; j<3; rp2[j]=R*u[j], ++j);
            --numPerturb;
            if (debug)
                mexPrintf("numPerturb=%d\n",numPerturb);
            brokenloop=false;
        }
        else if (brokenloop && numPerturb<=0) {
            st=255;
            forever=false;
        }   
    }
    
	return st;
    
    // Exit
}
