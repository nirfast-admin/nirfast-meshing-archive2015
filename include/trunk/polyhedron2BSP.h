/*
 *  polyhedron2BSP.h
 *  polyhedron2BSP
 *
 *  Created by Hamid_R_Ghadyani on 6/4/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __polyhedron2BSP
#define __polyhedron2BSP

#include "init.h"
#include "CPolygon.h"
#include "BSPNode.h"
#include "FileOperation.h"

#define PrinciplePlaneDelta 1.0

class Polyhedron2BSP {
	
public:
	Polyhedron2BSP();
	Polyhedron2BSP(std::vector<Polygon *> &inputpoly);
	Polyhedron2BSP(double *p, unsigned long *ele, unsigned long np, unsigned long ne, int nnpe);

	static ULONG maxdepth;
	std::vector<bool> polygonmarker;

	~Polyhedron2BSP();

	int ReadPolyhedronFromFile(std::string infilename);
	void InitFromMatlabMex(double *p, unsigned long *ele, unsigned long np, unsigned long ne, int nnpe);

	BSPNode* GetBSP_SolidLeaf_no_split();
	int IsInside(Point& p, double PlaneTHK);
	void SetBBX(double BBX[6]);
	void SetInputPolyhedron(std::vector<Polygon *> &inputpoly);
	//void SetInputPolyhedron(std::vector<Point *> verts, std::vector<unsigned long>

	double GetPlaneThickness();
	int DeleteTree();

	std::vector<Point *> _points;
	std::vector<unsigned long> _ele;
private:
	BSPNode* _BuildBSPTree_SL_NS(std::vector<Polygon *> &polygons, unsigned long depth, int label);
	Plane3D PickSplittingPlane(std::vector<Polygon *> &polygons, unsigned long depth);
	int PointInSolidSpace(BSPNode *node, Point& p, double PlaneTHK = TinyZero);
	BSPNode* _delete_node(BSPNode *);
	BSPNode* _root;
	std::vector<Polygon *> _inputpoly;
	std::vector<Plane3D> RequiredSplitPlanes;
	void SetupRequiredSplitPlanes();
	double* _bbx; //[xmin ymin zmin xmax ymax zmax]
	double _mindist;
	bool _isbuilt;
	REAL macheps;
	
};

extern REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
extern REAL exactinit();
#endif
