#include "polyhedron2BSP.h"

ULONG Polyhedron2BSP::maxdepth = 0;

Polyhedron2BSP::Polyhedron2BSP() : _isbuilt(false) {
	this->macheps = exactinit();
}
Polyhedron2BSP::Polyhedron2BSP(std::vector<Polygon *> &inputpoly) {
	this->macheps = exactinit();
	this->_inputpoly = inputpoly;
	this->polygonmarker.assign(this->_inputpoly.size(),false);
	SetupRequiredSplitPlanes();
}
Polyhedron2BSP::Polyhedron2BSP(double *p, unsigned long *ele, unsigned long np, unsigned long ne, int nnpe) {
	this->macheps = exactinit();
	this->_inputpoly.clear();

	for (unsigned long i=0; i<ne; ++i) {
		std::vector<Point *> fooverts;
		for (int j=0; j<nnpe; ++j) {
			unsigned long nodeid = (ele[i*nnpe+j]-1) * 3; // *3 for three coordinate components
			fooverts.push_back(new Point(p[nodeid], p[nodeid+1], p[nodeid+2]));
		}
		this->_inputpoly.push_back(new Polygon(fooverts));
	}
	_bbx = new double[6];
	for (unsigned long i=0; i<np; ++i) {
		if (i==0) {
			for (int j=0; j<6; j=j+3) {
				this->_bbx[j]   = p[i*3];
				this->_bbx[j+1] = p[i*3+1];
				this->_bbx[j+2] = p[i*3+2];
			}
		}
		if (p[i*3] < _bbx[0]) _bbx[0] = p[i*3];
		if (p[i*3] > _bbx[3]) _bbx[3] = p[i*3];

		if (p[i*3+1] < _bbx[1]) _bbx[1] = p[i*3+1];
		if (p[i*3+1] > _bbx[4]) _bbx[4] = p[i*3+1];

		if (p[i*3+2] < _bbx[2]) _bbx[2] = p[i*3+2];
		if (p[i*3+2] > _bbx[5]) _bbx[5] = p[i*3+2];
	}
	double x,y,z,longest;
	x = _bbx[3] - _bbx[0];
	y = _bbx[4] - _bbx[1];
	z = _bbx[5] - _bbx[2];
	longest = sqrt(x * x + y * y + z * z);
	if (longest == 0.0) {
		std::cerr << "\n\tError:  The point set is trivial.\n";
	    exit(1);
	}
	// Two identical points are distinguished by '_mindist'.
	_mindist = longest * this->macheps;

	// Setup initial split planes
	this->polygonmarker.assign(this->_inputpoly.size(),false);
	SetupRequiredSplitPlanes();
}

// This is called from a Matlab mex file. Since Matlab stores variables in a column based format we need this function.
void Polyhedron2BSP::InitFromMatlabMex(double *p, unsigned long *ele, unsigned long np, unsigned long ne, int nnpe) {
	this->macheps = exactinit();
	this->_inputpoly.clear();

	for (unsigned long i=0; i<ne; ++i) {
		std::vector<Point *> fooverts;
		for (int j=0; j<nnpe; ++j) {
			unsigned long nodeid = (ele[i + ne*j]-1);
			fooverts.push_back(new Point(p[nodeid], p[nodeid + np], p[nodeid + np + np]));
			this->_points.push_back(fooverts.back());
			this->_ele.push_back(nodeid);
		}
		Polygon *foop = new Polygon(fooverts);
		foop->id = i+1;
		this->_inputpoly.push_back(foop);
	}
	
	_bbx = new double[6];
	for (unsigned long i=0; i<np; ++i) {
		if (i==0) {
			for (int j=0; j<6; j=j+3) {
				this->_bbx[j]   = p[i]; // p[i + 0*np]
				this->_bbx[j+1] = p[i + np];
				this->_bbx[j+2] = p[i + np + np];
			}
		}
		if (p[i] < _bbx[0]) _bbx[0] = p[i]; // p[i + 0*np]
		if (p[i] > _bbx[3]) _bbx[3] = p[i]; // p[i + 0*np]

		if (p[i + np] < _bbx[1]) _bbx[1] = p[i + np];
		if (p[i + np] > _bbx[4]) _bbx[4] = p[i + np];

		if (p[i + np + np] < _bbx[2]) _bbx[2] = p[i + np + np];
		if (p[i + np + np] > _bbx[5]) _bbx[5] = p[i + np + np];
	}
	double x,y,z,longest;
	x = _bbx[3] - _bbx[0];
	y = _bbx[4] - _bbx[1];
	z = _bbx[5] - _bbx[2];
	longest = sqrt(x * x + y * y + z * z);
	if (longest == 0.0) {
		std::cerr << "\n\tError:  The point set is trivial.\n";
	    exit(1);
	}
	// Two identical points are distinguished by '_mindist'.
	_mindist = longest * this->macheps;

	_inputpoly[0]->GetPlane()->plane_thk_epsilon = _mindist;
	// Setup initial split planes
	SetupRequiredSplitPlanes();
	this->polygonmarker.assign(this->_inputpoly.size(),false);
}

void Polyhedron2BSP::SetInputPolyhedron(std::vector<Polygon *> &inputpoly) {
	
	this->_inputpoly = inputpoly;
	this->polygonmarker.assign(this->_inputpoly.size(),false);
}

void Polyhedron2BSP::SetBBX(double BBX[]) {
	double x,y,z,longest;
	_bbx = new double[6]; 
	for (int i=0; i<6; _bbx[i] = BBX[i], ++i);

	x = _bbx[3] - _bbx[0];
	y = _bbx[4] - _bbx[1];
	z = _bbx[5] - _bbx[2];
	longest = sqrt(x * x + y * y + z * z);
	if (longest == 0.0) {
		std::cerr << "\n\tError:  The point set is trivial.\n";
	    exit(1);
	}
	// Two identical points are distinguished by '_mindist'.
	_mindist = longest * this->macheps;
	SetupRequiredSplitPlanes();
}
void Polyhedron2BSP::SetupRequiredSplitPlanes() {
	// Most brain/breast shapes can be considered convex, so we also add 3 more planes parallel to XY/XZ/YZ planes that pass through centroid of object
	Point A, B;
	Point cen((_bbx[0]+_bbx[3])/2., (_bbx[1]+_bbx[4])/2., (_bbx[2]+_bbx[5])/2.);

	

	Plane3D foo1;
	double foomindist = 0.1;
	// We require 6 BBX planes to be part of split plane set
	Point llc(_bbx[0] - foomindist, _bbx[1] - foomindist, _bbx[2] - foomindist);
	Point urc(_bbx[3] + foomindist, _bbx[4] + foomindist, _bbx[5] + foomindist);

	
	A = urc; B = urc;
	A.y = urc.y - PrinciplePlaneDelta; B.z = urc.z - PrinciplePlaneDelta;
	foo1 = Plane3D(A,B,urc);
	this->RequiredSplitPlanes.push_back(foo1);

	A = urc; B = urc;
	A.z = urc.z - PrinciplePlaneDelta; B.x = urc.x - PrinciplePlaneDelta;
	foo1 = Plane3D(A,B,urc);
	this->RequiredSplitPlanes.push_back(foo1);

	A = llc; B = llc;
	A.z = llc.z + PrinciplePlaneDelta; B.y = llc.y + PrinciplePlaneDelta;
	foo1 = Plane3D(A,B,llc);
	this->RequiredSplitPlanes.push_back(foo1);

	A = llc; B = llc; 
	A.x = llc.x + PrinciplePlaneDelta; B.z = llc.z + PrinciplePlaneDelta;
	foo1 = Plane3D(A,B,llc);
	this->RequiredSplitPlanes.push_back(foo1);

	A = urc; B = urc; 
	A.x = urc.x - PrinciplePlaneDelta; B.y = urc.y - PrinciplePlaneDelta;
	foo1 = Plane3D(A,B,urc);
	this->RequiredSplitPlanes.push_back(foo1);

	A = llc; B = llc; 
	A.y = llc.y + PrinciplePlaneDelta; B.x = llc.x + PrinciplePlaneDelta;
	foo1 = Plane3D(A,B,llc);
	this->RequiredSplitPlanes.push_back(foo1);

	/////////////////////////////////////////////
	A.x = cen.x; A.y = cen.y + PrinciplePlaneDelta; A.z = cen.z;
	B.x = cen.x - PrinciplePlaneDelta; B.y = cen.y; B.z = cen.z;
	foo1 = Plane3D(A,B,cen);
	this->RequiredSplitPlanes.push_back(foo1);

	A.x = cen.x; A.y = cen.y; A.z = cen.z - PrinciplePlaneDelta;
	B.x = cen.x - PrinciplePlaneDelta; B.y = cen.y; B.z = cen.z;
	foo1 = Plane3D(A,B,cen);
	this->RequiredSplitPlanes.push_back(foo1);

	A.x = cen.x; A.y = cen.y; A.z = cen.z - PrinciplePlaneDelta;
	B.x = cen.x; B.y = cen.y + PrinciplePlaneDelta; B.z = cen.z;
	foo1 = Plane3D(A,B,cen);
	this->RequiredSplitPlanes.push_back(foo1);
	/*this->RequiredSplitPlanes.push_back(foo1);
	temp.n = Vector(0.,0.,-1.); // min. XY plane
	temp.d = -(_bbx[2] - _mindist);
	this->RequiredSplitPlanes.push_back(temp);
	temp.n = Vector(0.,0.,1.); // max. XY plane
	temp.d = -(_bbx[5] + _mindist);
	this->RequiredSplitPlanes.push_back(temp);
	temp.n = Vector(0.,-1.,0.); // min. XZ plane
	temp.d = -(_bbx[1] - _mindist);
	this->RequiredSplitPlanes.push_back(temp);
	temp.n = Vector(0.,1.,0.); // max. XZ plane
	temp.d = -(_bbx[4] + _mindist);
	this->RequiredSplitPlanes.push_back(temp);
	temp.n = Vector(-1.,0.,0.); // min. YZ plane
	temp.d = -(_bbx[0] - _mindist);
	this->RequiredSplitPlanes.push_back(temp);
	temp.n = Vector(1.,0.,0.); // max. YZ plane
	temp.d = -(_bbx[3] + _mindist);
	this->RequiredSplitPlanes.push_back(temp);*/
		
}

Polyhedron2BSP::~Polyhedron2BSP() {
	delete[] _bbx;
	for (ULONG i=0; i<_points.size(); delete _points[i], ++i);
	_inputpoly.clear();
}


BSPNode* Polyhedron2BSP::GetBSP_SolidLeaf_no_split() {
	if (_isbuilt)
		return this->_root;
	else {
		if (this->_inputpoly.empty())
			return (BSPNode *) NULL;
		srand( (unsigned)time( NULL ) ); // Initialize the seed for random generator.
		_root = _BuildBSPTree_SL_NS(this->_inputpoly, 0, BSPNode::OUT);
		_isbuilt = true;
		return _root;
	}
}

BSPNode* Polyhedron2BSP::_BuildBSPTree_SL_NS(std::vector<Polygon *> &polygons, unsigned long depth, int label) {

    // Get number of polygons in the input vector
    ULONG numPolygons = polygons.size();
	Plane3D splitPlane;

	if (polygons.empty()) {
		return new BSPNode(label);
	}

	std::vector<Polygon *> frontList, backList;

	// Using following section, splitting planes are chosen based on a scoring system.
	///////////////////////////////////////////////////////////////
	if (depth > maxdepth) maxdepth = depth;
    // Select best possible partitioning plane based on the input geometry
	// if (numPolygons != 1)
	splitPlane = PickSplittingPlane(polygons, depth);
	/*else {
		splitPlane = Plane3D(polygons[0]->GetVertexPtr(0), polygons[0]->GetVertexPtr(1), polygons[0]->GetVertexPtr(2));
		// if (this->polygonmarker[ polygons[0]->id - 1])
		// 	std::cout << "  _BuildBSPTree: polygon's plane has already been used!" << std::endl;
		this->polygonmarker[ polygons[0]->id - 1] = true;
		BSPNode *frontTree = _BuildBSPTree_SL_NS(frontList, depth+1, BSPNode::OUT);
		BSPNode *backTree = _BuildBSPTree_SL_NS(backList, depth+1, BSPNode::IN);
		return new BSPNode(frontTree, backTree, splitPlane, depth);
	}*/

    // Test each polygon against the dividing plane, adding them
    // to the front list, back list, or both, as appropriate
    for (unsigned long i = 0; i < numPolygons; i++) {
        Polygon *poly = polygons[i]; //*frontPart, *backPart;
		//int st;
		switch (poly->ClassifyPolygonToPlane(splitPlane)) {
			case Polygon::POLYGON_COPLANAR_WITH_PLANE:
				// What's done in this case depends on what type of tree is being
				// built. For a node-storing tree, the polygon is stored inside
				// the node at this level (along with all other polygons coplanar
				// with the plane). Here, for a leaf-storing tree, coplanar polygons
				// are sent to either side of the plane. In this case, to the front
				// side, by falling through to the next case
				// frontList.push_back(poly);
				//if (numPolygons==1)
				this->polygonmarker[poly->id - 1] = true;
				break;
			case Polygon::POLYGON_IN_FRONT_OF_PLANE:
				frontList.push_back(poly);
				break;
			case Polygon::POLYGON_BEHIND_PLANE:
				backList.push_back(poly);
				break;
			case Polygon::POLYGON_STRADDLING_PLANE:
				// Split polygon to plane and send a part to each side of the plane
				//SplitPolygon(*poly, splitPlane, &frontPart, &backPart);
				frontList.push_back(poly);
				backList.push_back(poly);
				break;
		}
	}

    // Recursively build child subtrees and return new tree root combining them
	BSPNode *frontTree = _BuildBSPTree_SL_NS(frontList, depth+1, BSPNode::OUT);
	BSPNode *backTree = _BuildBSPTree_SL_NS(backList, depth+1, BSPNode::IN);
    return new BSPNode(frontTree, backTree, splitPlane, depth);
}

Plane3D Polyhedron2BSP::PickSplittingPlane(std::vector<Polygon *> &polygons, unsigned long depth)
{
    // Blend factor for optimizing for balance or splits (should be tweaked)
    const float K = 0.8f;
    // Variables for tracking best splitting plane seen so far
    Plane3D bestPlane;
	// Variable for tracking index of the best plane
	ULONG idx = 0;

	///////////////////////////////////////////////////////////////////////////////////
	/*if (this->polygonmarker[ polygons[0]->id - 1])
			std::cout << "  PickSplittingPlane: polygon's plane has already been used!" << std::endl;*/
	srand( (unsigned)time( NULL ) );
	idx = myrand((ULONG) polygons.size());
	assert(idx<(ULONG)polygons.size() && idx>=0);
	bestPlane = *(polygons[idx]->GetPlane());
	this->polygonmarker[ polygons[idx]->id - 1] = true;
	return bestPlane;
	///////////////////////////////////////////////////////////////////////////////////

	float bestScore = std::numeric_limits<float>::max();

	/*if (depth < this->RequiredSplitPlanes.size()) { // We should first split using required splitting planes
		return this->RequiredSplitPlanes[depth];
	}*/
    // Try the plane of each polygon as a dividing plane
    for (unsigned long i = 0; i < polygons.size(); i++) {
        unsigned long numInFront = 0, numBehind = 0, numStraddling = 0;
        Plane3D* plane = polygons[i]->GetPlane();
        // Test against all other polygons
        for (unsigned long j = 0; j < polygons.size(); j++) {
            // Ignore testing against self
            if (i == j) continue;
            // Keep standing count of the various poly-plane relationships
			switch (polygons[j]->ClassifyPolygonToPlane(plane)) {
			case Polygon::POLYGON_COPLANAR_WITH_PLANE:
                /* Coplanar polygons treated as being in front of plane */
            case Polygon::POLYGON_IN_FRONT_OF_PLANE:
                numInFront++;
                break;
            case Polygon::POLYGON_BEHIND_PLANE:
                numBehind++;
                break;
            case Polygon::POLYGON_STRADDLING_PLANE:
                numStraddling++;
                break;
            }
        }
        // Compute score as a weighted combination (based on K, with K in range
        // 0..1) between balance and splits (lower score is better)
		float score = K * numStraddling + (1.0f - K) * std::fabs((float)((float)numInFront - (float)numBehind));
        if (score < bestScore) {
            bestScore = score;
            bestPlane = *plane;
			idx = i;
        }
    }
	/*if (this->polygonmarker[ polygons[0]->id - 1])
			std::cout << "  PickSplittingPlane: polygon's plane has already been used!" << std::endl;*/
	this->polygonmarker[ polygons[idx]->id - 1] = true;
    return bestPlane;
}

int Polyhedron2BSP::IsInside(Point& p, double PlaneTHK) {
	
	BSPNode *node = this->GetBSP_SolidLeaf_no_split();
	return this->PointInSolidSpace(node, p, PlaneTHK);
}

int Polyhedron2BSP::PointInSolidSpace(BSPNode *node, Point& p, double PlaneTHK)
// Using a solid-leaf BSP, determines if point p is inside, outside
// or on the boundary of polyhedron
// Returns:
// 0 : Outside
// 1 : Inside
// 2 : On the boundary
{
    while (!node->IsLeaf()) {
		if (true) {
            double a[3],b[3],c[3],d[3];
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

            if (ret < 0.0)
                node = node->frontnode;
            else if (ret > 0.0)
                node = node->backnode;
            else if (ret == 0.0) {
                int front = PointInSolidSpace(node->frontnode, p);
                int back  = PointInSolidSpace(node->backnode,  p);
                // If results agree, return that, else point is on boundary
                //std::cout << " Query point on node! " << ((front == back) ? front : 2) << std::endl;
                //std::cout << " P = (" << p.x << ", " << p.y << ", " << p.z << ')' << std::endl;
                return (front == back) ? front : 2;
            }
            else {
                std::cout << std::endl << 
                "  There seems to be bug!!!" << std::endl << std::endl;
            }
        }
        else {
            // Compute distance of point to dividing plane
            double dist = node->myplane.n*p + node->myplane._d;
            if (dist > PlaneTHK) {
                // Point in front of plane, so traverse front of tree
                node = node->frontnode;
            } else if (dist < -PlaneTHK) {
                // Point behind of plane, so traverse back of tree
                node = node->backnode;
            } else {
                // Point on dividing plane; must traverse both sides
                int front = PointInSolidSpace(node->frontnode, p);
                int back  = PointInSolidSpace(node->backnode,  p);
                // If results agree, return that, else point is on boundary
                return (front == back) ? front : 2;
            }
        }
             
             
    }
    // Now at a leaf, inside/outside status determined by solid flag
    return node->IsSolid() ? 1 : 0;    
}

int Polyhedron2BSP::ReadPolyhedronFromFile(std::string infn) {
	CFileOperation myfile;
	std::ifstream ifs;
	std::string nodefn, elefn;

	this->_inputpoly.clear();
	// Initiate exact arithmetic
	this->macheps = exactinit();
	std::string ext = myfile.GetExtension(infn);
	size_t elefound, nodefound;
	elefound =  ext.rfind("ele"); nodefound = ext.rfind("node");
	if (elefound==std::string::npos && nodefound==std::string::npos) { // need ele/node extensions
		std::cerr << std::endl << "  Polyhedron2BSP::ReadPolyhedronFromFile: Need .ele/.node filenames!" << std::endl;
		return 1;
	}
	if (elefound==std::string::npos) {
		nodefn = infn;
		elefn = myfile.MakeFilename(infn, ".node", ".ele");
	}
	else {
		nodefn = myfile.MakeFilename(infn, ".ele", ".node");
		elefn = infn;
	}
	MeshIO polyio(nodefn, elefn);
	/*std::vector<Point *> points;
	std::vector<unsigned long> ele; // We assume that we are going to read a triangular faceted polyhedron
	*/
	int nnpe;
	if ((nnpe = polyio.ReadPolyhedron(_points, _ele)) < 0) {
		std::cerr << std::endl << "\tError in reading node/ele file!" << std::endl;
		return 2;
	}

	// Populate the input polyhedron
	for (unsigned long i=0; i< _ele.size()/nnpe; ++i) {
		std::vector<Point *> fooverts;
		// get points of each polygonal face of the polyhedron
		for (int j=0; j<nnpe; ++j) {
			fooverts.push_back(_points[ _ele[i*nnpe+j] - 1 ]);
		}
		Polygon *foop = new Polygon(fooverts);
		foop->id = i+1;
		this->_inputpoly.push_back(foop);
	}
	// Get bbx of the polyhedron
	_bbx = new double[6];
	for (unsigned long i=0; i<_points.size(); ++i) {
		if (i==0) {
			for (int j=0; j<6; j=j+3) {
				this->_bbx[j]   = _points[i]->x;
				this->_bbx[j+1] = _points[i]->y;
				this->_bbx[j+2] = _points[i]->z;
			}
		}
		if (_points[i]->x < _bbx[0]) _bbx[0] = _points[i]->x;
		if (_points[i]->x > _bbx[3]) _bbx[3] = _points[i]->x;

		if (_points[i]->y < _bbx[1]) _bbx[1] = _points[i]->y;
		if (_points[i]->y > _bbx[4]) _bbx[4] = _points[i]->y;

		if (_points[i]->z < _bbx[2]) _bbx[2] = _points[i]->z;
		if (_points[i]->z > _bbx[5]) _bbx[5] = _points[i]->z;
	}
	double x,y,z,longest;
	x = _bbx[3] - _bbx[0];
	y = _bbx[4] - _bbx[1];
	z = _bbx[5] - _bbx[2];
	longest = sqrt(x * x + y * y + z * z);
	if (longest == 0.0) {
		std::cerr << "\n\tError:  The point set is trivial.\n";
	    exit(1);
	}
	// Two identical points are distinguished by '_mindist'.
	_mindist = longest * this->macheps;

	_inputpoly[0]->GetPlane()->plane_thk_epsilon = this->_mindist;
	// Setup initial split planes
	SetupRequiredSplitPlanes();
	this->polygonmarker.assign(this->_inputpoly.size(),false);
	return 0;
}

int Polyhedron2BSP::DeleteTree() {
	if (!_isbuilt) {
		return 0;
	}
	else {
		_delete_node(this->_root);
		delete _root;
		_isbuilt = false;
		return 0;
	}
}
BSPNode* Polyhedron2BSP::_delete_node(BSPNode* node) {
	if (node->IsLeaf()) {
		return node;
	}
	else {
		BSPNode *ret = _delete_node(node->frontnode);
		delete ret;
		ret = _delete_node(node->backnode);
		delete ret;
		return node;
	}
}

double Polyhedron2BSP::GetPlaneThickness() {
	return this->_inputpoly[0]->GetPlane()->plane_thk_epsilon;
}
