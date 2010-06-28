#include "polyhedron2BSP.h"
#include "FileOperation.h"
#include <iostream>
#include <fstream>
#include <time.h>

#include "mex.h"

int PointInSolidSpace(BSPNode *node, Point& p, double PlaneTHK = TinyZero);
std::vector< std::vector<int> > core(std::vector<std::string>&surface_fnames);

/*
This program detects the the inclusion of nested 3D surfaces and how they are nested.
This is required to formulate BEM matrices.

This MEX file should be called from Matlab using following format:

>> surface_relations_mex(fnprefix,no_regions);

'fnprefix' is first part of .node/.ele filenames: 'fnprefix_4.node/.ele'
'no_regions' is total number of sepearate homogenous regions.

Output: a text file "surface_relations.txt". For exmaple if we have 5 regions, this text file will be:

1 2 3 4
2 0
3 5
4 0
5 0

meaning that surface 1 has three direct children (2, 3 and 4). surface 3 has only one child (5) and
the rest of surfaces don't have any child.

If any pair of surfaces are intersecting, the program will not produce "surface_relations.txt" and
will return a value of '5'. Otherwise it will return 0.

*/

/* How to compile:

Windows:
mex -v -DCPU86 -DWIN32 -I./meshlib surface_relations_mex.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp

Linux:
mex -v -DCPU86 -DLinux -I./meshlib ssurface_relations_mex.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp

Mac OSX:
mex -v -I./meshlib                 surface_relations_mex.cpp meshlib/polyhedron2BSP.cpp meshlib/CPoint.cpp meshlib/CVector.cpp meshlib/Plane3D.cpp meshlib/BSPNode.cpp meshlib/MeshIO.cpp meshlib/FileOperation.cpp meshlib/CPolygon.cpp meshlib/predicates.cpp meshlib/CStopWatch.cpp
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs!=2)
		mexErrMsgTxt("surface_relations_mex needs 2 input arguments!");
	/*if (nlhs!=1)
		mexErrMsgTxt("surface_relations_mex needs 1 output variable!");*/
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("first input argument to surface_relations_mex should filename prefix of .inp files!");

	char *fname_prefix = new char[256];
    // Get filename prefix for all the .inp files
	fname_prefix = mxArrayToString(prhs[0]);
    // Get number of regions
	int no_regions = (int)(*mxGetPr(prhs[1]));

	std::vector<std::string> snames;
	char buffer[256];
	std::string tmpstring("");
    
	for (int i=0; i<no_regions; ++i) {
#if defined(WIN32) || defined(WIN64) || defined(WINDOWS) || defined(windows)
		tmpstring = (((std::string)fname_prefix) + *(_itoa(i+1,buffer,10)));
#else
        sprintf(buffer,"%d",i+1);
        tmpstring = ((std::string)fname_prefix) + *buffer;
#endif
		tmpstring = tmpstring + ".node";
		snames.push_back(tmpstring);
	}
	std::vector< std::vector<int> > relations = core(snames);
}

std::vector< std::vector<int> > core(std::vector<std::string>&surface_fnames) {
	int c = 0;
	CFileOperation myfile;

	int NumberOfDomains = (int) surface_fnames.size();
	Polyhedron2BSP *tree = new Polyhedron2BSP[NumberOfDomains];
	std::vector< std::vector<int> > relations;

	for (int i=0; i<NumberOfDomains; ++i) {
		int j = tree[i].ReadPolyhedronFromFile(surface_fnames[i]);
		if (j!=0) {
			mexPrintf("surface_relations_mex: Could not read one of the input files!\n");
//			delete [] tree;
			return relations;
		}
		tree[i].SetSplitType(2); // Full scoring system
	}


	std::vector<BSPNode *> root;
	root.resize(NumberOfDomains);
	mexPrintf("Constructiong the BSP tree... ");
	double tic = CPU_TIME;
	for (int i=0; i<NumberOfDomains; ++i) {
		root[i] = tree[i].GetBSP_SolidLeaf_no_split();
	}
	double toc = CPU_TIME;
	mexPrintf("done (time: %.12g).\n",(toc-tic)/(double)CLOCKS_PER_SEC);

	std::ofstream ofs;
	std::string outfn("surface_relations.txt");
	std::vector<int> foo1, foo2;

	int mymax;
	std::vector<int> rank(NumberOfDomains,0);
	std::vector<bool> idxmarker(NumberOfDomains,false);
	std::vector< std::vector<int> >::iterator ite;

	bool allin_flag = true, inflag = false, loopflag = true;;
	
	relations.assign(NumberOfDomains,0);
	std::vector<int> foo;
	for (int i=0; i<NumberOfDomains; ++i) {
		double foothk = tree[i].GetPlaneThickness();
		for (int j=0; j<NumberOfDomains; ++j) {
			if (i==j) continue;
			for (ULONG ii=0; ii<tree[j]._points.size(); ++ii) {
				Point P = *(tree[j]._points[ii]);
				//int ret = PointInSolidSpace(root[i], P, foothk);
				int ret = tree[i].IsInside(P, foothk);
				if (ret == 0) {
					allin_flag = false;
					// break;
				}
				else if (ret == 1 || ret == 2) {
					inflag = true;
				}
			}
			if (allin_flag == false && inflag == true) { // intersection surfaces, not suitable for BEM
				mexPrintf("  Warning!\n");
				mexPrintf("   This set of surfaces can not be used in BEM formulation!\n");
				mexPrintf("   Reason: they are NOT mutually exclusive sub-surfaces!\n\n");
				allin_flag = true; inflag = false;
				goto exit;
			}
			else if (allin_flag == true) {
				//std::cout << "just before bug" << std::endl;
				foo = relations[i];
				foo.push_back(j);
				relations[i] = foo;
				inflag = false;
			}
			else if (allin_flag == false && inflag == false) {
				allin_flag = true;
			}
		}
	}
	// Sort the relations vector based on number of childern (max to min)
	c=0;
	while (loopflag) {
		loopflag = false;
		for (int i=0; i<relations.size(); ++i) {
			if (idxmarker[i]==true)
				continue;
			else {
				mymax = i;
				loopflag = true;
			}
			for (int j=i+1; j<relations.size(); ++j) {
				if (idxmarker[j]==true)
					continue;
				if (relations[j].size() > relations[mymax].size())
					mymax = j;
			}
			idxmarker[mymax]=true;
			rank[c] = mymax;
			++c;
		}
	}

	for (int i=0; i<relations.size()-1; ++i) {
		foo1 = relations[rank[i]];
		for (int j=i+1; j<relations.size(); ++j) {
			foo2 = relations[rank[j]];
			for (int ii=0; ii<foo2.size(); ++ii) {
				for (std::vector<int>::iterator jj=foo1.begin(); jj!=foo1.end(); ++jj) {
					if (foo2[ii]==*jj) {
						foo1.erase(jj);
						relations[rank[i]] = foo1;
						break;
					}
				}
			}
		}
	}
	
	myfile.OpenOutFile(ofs,outfn,"text");
	c=1;
	// Write a text file explaining the relationships between surfaces
	for (int i=0; i<NumberOfDomains; ++i) {
		foo1 = relations[rank[i]];
		ofs << rank[i] + 1;
		if (!foo1.empty()) {
			for (int j=0; j<foo1.size(); ++j) {
				ofs << ' ' << foo1[j] + 1;
			}
			ofs << std::endl;
		}
		else {
			ofs << " 0" << std::endl;
		}
	}
	ofs.close();

exit:
	for (int i=0; i<NumberOfDomains; ++i) {
		tree[i].DeleteTree();
	}
	root.clear();
	return relations;
}

int PointInSolidSpace(BSPNode *node, Point& p, double PlaneTHK)
// Using a solid-leaf BSP, determines if point p is inside, outside
// or on the boundary of polyhedron defined by 'node'
// Returns:
// 0 : Outside
// 1 : Inside
// 2 : On the boundary
{
    while (!node->IsLeaf()) {
        // Compute distance of point to dividing plane
        //float dist = Dot(node->plane.n, p) - node->plane.d;
		double dist = node->myplane.n*p + node->myplane._d;

		/*double a[3],b[3],c[3],d[3];
		a[0] = node->myplane.GetThreePoints(0)->x;
		a[1] = node->myplane.GetThreePoints(0)->y;
		a[2] = node->myplane.GetThreePoints(0)->z;
		b[0] = node->myplane.GetThreePoints(1)->x;
		b[1] = node->myplane.GetThreePoints(1)->y;
		b[2] = node->myplane.GetThreePoints(1)->z;
		c[0] = node->myplane.GetThreePoints(2)->x;
		c[1] = node->myplane.GetThreePoints(2)->y;
		c[2] = node->myplane.GetThreePoints(2)->z;
		d[0] = p.x; d[1] = p.y; d[2] = p.z;
		double ret = orient3d(a,b,c,d);
		if (ret < -node->myplane.plane_thk_epsilon)
			node = node->frontnode;
		else if (ret > node->myplane.plane_thk_epsilon)
			node = node->backnode;
		else {
			int front = PointInSolidSpace(node->frontnode, p);
			int back = PointInSolidSpace(node->backnode, p);
            // If results agree, return that, else point is on boundary
            return (front == back) ? front : 2;
		}
	}
	return node->IsSolid() ? 1 : 0;*/
        if (dist > PlaneTHK) {
            // Point in front of plane, so traverse front of tree
			node = node->frontnode;
        } else if (dist < -PlaneTHK) {
            // Point behind of plane, so traverse back of tree
			node = node->backnode;
        } else {
            // Point on dividing plane; must traverse both sides
			int front = PointInSolidSpace(node->frontnode, p);
			int back = PointInSolidSpace(node->backnode, p);
            // If results agree, return that, else point is on boundary
            return (front == back) ? front : 2;
        }
    }
    // Now at a leaf, inside/outside status determined by solid flag
    return node->IsSolid() ? 1 : 0;
}
 
