#include "SurfaceParameterization.hpp"

JSurfaceParameterization :: JSurfaceParameterization()
{
    boundarytype = 0;
    weighttype   = 0;
    numIters     = 2000;
    gammaP       = 1.0;
    smooth       = 1;
    intrinsiclambda=0.5;
}
///////////////////////////////////////////////////////////////////////////////
int JSurfaceParameterization :: setPatchMesh( const JMeshPtr &m)
{
    patchMesh.reset();
    int val = m->getTopology()->isDisk();

    if( val == 0) {
        cout << "Warning: the patch is not topological disk " << endl;
        return 1;
    }
    patchMesh = m;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterization :: normalize( const JMeshPtr &mesh2d)
{
    if( mesh2d == nullptr) return;

    JMeshAffineTransform affine;
    affine.setMesh( mesh2d );
    affine.toCenter();
    affine.normalize();
}
///////////////////////////////////////////////////////////////////////////////
int JSurfaceParameterization :: writePly2File()
{
    if( patchMesh == nullptr) return 1;

    size_t numnodes = patchMesh->getSize(0);
    size_t numfaces = patchMesh->getSize(2);

    int numTriangles = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = patchMesh->getFaceAt(i);
        if( face->isActive() ) {
            int nn = face->getSize(0);
            if( nn == 3) numTriangles += 1;
            if( nn == 4) numTriangles += 2;
        }
    }

    ofstream ofile( "in.ply2", ios::out);
    if( ofile.fail()) return 2;

    ofile << numnodes  << endl;
    ofile << numTriangles << endl;

    map<JNodePtr,size_t> localID;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = patchMesh->getNodeAt(i);
        const Point3D &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
        localID[vtx] = i;
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = patchMesh->getFaceAt(i);
        if( face->isActive() ) {
            int nn = face->getSize(0);
            if( nn  == 3) {
                size_t  id0 = localID[face->getNodeAt(0)];
                size_t  id1 = localID[face->getNodeAt(1)];
                size_t  id2 = localID[face->getNodeAt(2)];
                ofile << "3 " << id0 << " " << id1 << " " << id2 << endl;
            }
            if( nn == 4) {
                size_t  id0 = localID[face->getNodeAt(0)];
                size_t  id1 = localID[face->getNodeAt(1)];
                size_t  id2 = localID[face->getNodeAt(2)];
                size_t  id3 = localID[face->getNodeAt(3)];
                ofile << "3 " << id0 << " " << id1 << " " << id2 << endl;
                ofile << "3 " << id0 << " " << id2 << " " << id3 << endl;
            }
        }
    }
    return 0;
}


///////////////////////////////////////////////////////////////////////////////

JMeshPtr JSurfaceParameterization :: readPly2File()
{
    ifstream ifile( "out.ply2", ios::in);
    if( ifile.fail()) return nullptr;

    JMeshPtr paramMesh = JMesh::newObject();

    size_t numnodes, numfaces;
    ifile >> numnodes;
    ifile >> numfaces;

    Point3D xyz;
    map<JNodePtr, JNodePtr> nodeMap;
    for( size_t i = 0; i < numnodes; i++) {
        ifile >> xyz[0] >> xyz[1] >> xyz[2];
        JNodePtr newnode = JNode::newObject();
        newnode->setXYZCoords(xyz);
        paramMesh->addObject(newnode);
        JNodePtr oldnode = patchMesh->getNodeAt(i);
        nodeMap[oldnode] = newnode;
    }

    numfaces = patchMesh->getSize(2);
    JNodeSequence nodes;

    JNodeSequence connect;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &srcface = patchMesh->getFaceAt(i);
        int nsize = srcface->getSize(0);
        connect.resize(nsize);
        for( int j = 0; j < nsize; j++)
            connect[j] = nodeMap[srcface->getNodeAt(j)];
        JFacePtr dstface = srcface->getClone();
        dstface->setNodes( connect );
        paramMesh->addObject(dstface);
    }
    normalize(paramMesh);
    return paramMesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JSurfaceParameterization :: LeastSquareConformalMap()
{
    JMeshEigenMatrix mat;

    cout << "Stage 1 " << endl;

    mat.setMesh(patchMesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    cout << "Stage 1 " << endl;

    // Fix two points on the boundary
    Eigen::VectorXi bnd, b(2,1);
    igl::boundary_loop(F, bnd);

    b(0) = bnd(0);
    b(1) = bnd(round(bnd.size()/2));

    Eigen::MatrixXd bc(2,2);
    bc<< 0,0,1,0;
    cout << "Stage 1 " << endl;

    Eigen::MatrixXd uv;
    // LSCM parametrization
    igl::lscm(V, F, b, bc, uv);
    cout << "Stage 1 " << endl;

    JMeshPtr paramMesh = mat.getMesh(uv, F);

    int val = 0;
    size_t numnodes = patchMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = patchMesh->getNodeAt(i);
        if( v0->hasAttribute("PartitionCorner") ) {
            const JNodePtr &v1 = paramMesh->getNodeAt(i);
            v1->setAttribute("PartitionCorner", val);
        }
    }
    cout << "Stage 1 " << endl;

    paramMesh->setName("Param_LSCM");
    normalize(paramMesh);
    return paramMesh;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JSurfaceParameterization :: AsRigidAsPossible()
{
    JMeshEigenMatrix mat;

    mat.setMesh(patchMesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    Eigen::VectorXi bnd;
    igl::boundary_loop(F,bnd);
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V,bnd,bnd_uv);

    Eigen::MatrixXd initial_guess;
    igl::harmonic(V,F,bnd,bnd_uv,1,initial_guess);

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 100;
    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V,F,2,b,arap_data);

    // Solve arap using the harmonic map as initial guess
    Eigen::MatrixXd uv = initial_guess;

    arap_solve(bc,arap_data,uv);

    JMeshPtr paramMesh = mat.getMesh(uv, F);

    int val = 0;
    size_t numnodes = patchMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = patchMesh->getNodeAt(i);
        if( v0->hasAttribute("PartitionCorner") ) {
            const JNodePtr &v1 = paramMesh->getNodeAt(i);
            v1->setAttribute("PartitionCorner",val);
        }
    }

    paramMesh->setName("Param_ARAP");
    normalize(paramMesh);
    return paramMesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JSurfaceParameterization :: getParamMesh()
{
    if( patchMesh == nullptr) return nullptr;

    if( weighttype == LEAST_SQUARE_CONFORMAL ) return LeastSquareConformalMap();
    if( weighttype == AS_RIGID_AS_POSSIBLE )   return AsRigidAsPossible();

    writePly2File();
    ostringstream oss;
    oss << "/home/csverma/Projects/QuadMesh/JaalMesh/external/Parameterization/surfparam ";
    oss << "-b" << boundarytype << " ";
    oss << "-g" << gammaP << " ";
    oss << "-l" << intrinsiclambda << " ";
    oss << "-n" << numIters << " ";
//  oss << "p " << paramtype;
    oss << "-s" << smooth << " ";
    oss << "-w" << weighttype << " ";
    oss << "-i " << "in.ply2" << " ";
    oss << "-o " << "out.ply2";

    int stat = system(oss.str().c_str() );
    if( stat == -1 ) return nullptr;
    JMeshPtr paramMesh = readPly2File();

    int val;
    size_t numnodes = patchMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = patchMesh->getNodeAt(i);
        if( v0->hasAttribute("PartitionCorner") ) {
            const JNodePtr &v1 = paramMesh->getNodeAt(i);
            v1->setAttribute("PartitionCorner",val);
        }
    }

    return paramMesh;
}

///////////////////////////////////////////////////////////////////////////////

