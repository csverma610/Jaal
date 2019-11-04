#include <fstream>
#include "Elasticity2D.hpp"
#include "MeshOptimization.hpp"
#include "MeshRefine.hpp"
#include "MeshDual.hpp"
#include "DelaunayMesh.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"

//extern int getTexCoords(const JFacePtr &face, const Point2D &xy, Point2D &uv);
//extern int isInside( const vector<Point2D> &polyPoints, const Point2D &qPoint);
//extern int getIntersection( const JFacePtr &f1, const JFacePtr &f2, vector<Point2D> &p);

int randomTangle, testmodel, checkTangle = 0, meshType = 4;
int elemOrder = 1, numGuassPnts, triGaussPnts;
int rotatemesh = 0, rotAngle = 0.0;

double youngModulus = 2.0E+11;
double poissonRatio = 0.25;

int  numGauss1    = 3;   // (1,3,7)
int  numGauss2    = 1;   // (1,3,7)
bool absJacobian  = 0;
int  fieldType    = 1;
int  shapeFamily  = 0;
int  checkConvex  = 1;

double circleCenter[2], circleRadius;
ofstream  outfile;

/////////////////////////////////////////////////////////////////////////////////////
void  SplitFace( const JMeshPtr &mesh, const JFacePtr &oldface)
{
   Point3D  xyz;
   Point2D uv;
   uv[0] = -0.8;
   uv[1] = -0.8;

   JNodeSequence v(5);
   QuadGeometry::getBilinearCoords(oldface, uv, xyz);
   v[4] = JVertex::newObject();
   v[4]->setXYZCoords(xyz);

   v[0] = oldface->getNodeAt(0);
   v[1] = oldface->getNodeAt(1);
   v[2] = oldface->getNodeAt(2);
   v[3] = oldface->getNodeAt(3);

   JFacePtr f0 = Quadrilateral::newObject( v[0], v[1], v[4], v[3] );
   JFacePtr f1 = Quadrilateral::newObject( v[1], v[2], v[3], v[4] );
   mesh->addObject(v[4] );
   mesh->addObject( f0 );
   mesh->addObject( f1 );
   oldface->setStatus( JMeshEntity::REMOVE);
}
/////////////////////////////////////////////////////////////////////////////////////

int SwapPositions( const JNodePtr &v1, const JNodePtr &v2, vector<size_t> &permute)
{
    const Point3D p1 = v1->getXYZCoords();
    const Point3D p2 = v2->getXYZCoords();

    v1->setXYZCoords(p2);
    v2->setXYZCoords(p1);

    int id1 = v1->getID();
    int id2 = v2->getID();
    permute[id1] = id2;
    permute[id2] = id1;
}

/////////////////////////////////////////////////////////////////////////////////////
void affineFunc( const Point3D &p, double &u, double &v)
{
    double x = p[0];
    double y = p[1];

    if( fieldType == 1 ) {
        u =  0.1*x - 0.2*y + 1.0E-10;
        v = -0.3*x + 0.5*y + 1.0E-10;
    } else {
        double a1  = 10;
        double a2  = 10;
        double nu = 0.3;
        u = a1*x*x - 4.0*a2*x*y/(1+nu);
        v = -4.0*a1*x*y/(1+nu) + a2*y*y;
    }
}
/////////////////////////////////////////////////////////////////////////////////////

double patchTest( const JMeshPtr &mesh)
{
    JElasticity2D  elastic;
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle(checkTangle);
    elastic.checkConvexity(0);
    elastic.useAbsoluteJacobian(absJacobian);
    elastic.setOrder(elemOrder);
    elastic.setNumFaceGaussPoints(numGauss1);
    elastic.setNumTangleFaceGaussPoints(numGauss2);
    elastic.setMesh(mesh);

/*
    int etype = mesh->getTopology()->isHomogeneous(2);
    JNodeSequence boundnodes;
    mesh->getTopology()->getBoundary(boundnodes);

    assert( !boundnodes.empty() ) ;

    double uexact, vexact;
    for( int i = 0; i < boundnodes.size(); i++) {
        const Point3D &xyz =  boundnodes[i]->getXYZCoords();
        affineFunc(xyz, uexact,vexact);
        elastic.fixX( boundnodes[i], uexact);
        elastic.fixY( boundnodes[i], vexact);
    }

    elastic.setShapeFamily(shapeFamily);
    elastic.solve();

    vector<double> usol, vsol;
    elastic.getDisplacements(usol, vsol);

    size_t numnodes = mesh->getSize(0);
    cout << fixed << setprecision(15);
    double maxu = 0.0;
    double maxv = 0.0;

    vector<double> error;
    error.reserve( numnodes);
    for( int i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( !vtx->isBoundary() )  {
        const Point3D &xyz =  mesh->getNodeAt(i)->getXYZCoords();
        affineFunc(xyz, uexact,vexact);
        double du = usol[i] - uexact;
        double dv = vsol[i] - vexact;
        error[i] = sqrt( du*du + dv*dv);
        if( fabs(uexact) > 1.0E-10)
            maxu = max(maxu, fabs(du)/(fabs(uexact)));
        else
            maxu = max(maxu, fabs(du));

        if( fabs(vexact) > 1.E-10)
            maxv = max(maxv, fabs(dv)/(fabs(vexact)));
        else
            maxv = max(maxv, fabs(dv));
        }
    }

    cout << scientific << endl;
    double maxerror = *max_element( error.begin(), error.end() );
    cout << "Maximum relative errors " << 100*maxu << "  " << 100*maxv << endl;
    outfile << max(100*maxu, 100*maxv) << endl;
    return maxerror;
*/
}
///////////////////////////////////////////////////////////////////////////////
#ifdef CSV

void testRotateQuadPatch1()
{
    JMeshPtr mesh = JMesh::newObject();

    JNodeSequence nodes(8);
    for( int i = 0; i < 8; i++) {
        nodes[i] = JVertex::newObject();
        nodes[i]->setID(i);
    }
    Point3D xyz;
    xyz[2] = 0.0;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    nodes[0]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    nodes[1]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    nodes[2]->setXYZCoords(xyz);

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    nodes[3]->setXYZCoords(xyz);

    xyz[0] = 1.0/3.0;
    xyz[1] = 1.0/3.0;
    nodes[4]->setXYZCoords(xyz);

    xyz[0] = 2.0/3.0;
    xyz[1] = 1.0/3.0;
    nodes[5]->setXYZCoords(xyz);

    xyz[0] = 2.0/3.0;
    xyz[1] = 2.0/3.0;
    nodes[6]->setXYZCoords(xyz);

    xyz[0] = 1.0/3.0;
    xyz[1] = 2.0/3.0;
    nodes[7]->setXYZCoords(xyz);

    mesh->addObjects(nodes);

    JFaceSequence faces(5);
    faces[0] = Quadrilateral::newObject( nodes[0], nodes[1], nodes[5], nodes[4] );
    faces[1] = Quadrilateral::newObject( nodes[1], nodes[2], nodes[6], nodes[5] );
    faces[2] = Quadrilateral::newObject( nodes[2], nodes[3], nodes[7], nodes[6] );
    faces[3] = Quadrilateral::newObject( nodes[3], nodes[0], nodes[4], nodes[7] );
    faces[4] = Quadrilateral::newObject( nodes[4], nodes[5], nodes[6], nodes[7] );
    mesh->addObjects( faces);

    if( meshType == 3)  {
        JMeshPtr trimesh = AllTriMeshGenerator::fromQuadMesh(mesh,4);
        mesh = trimesh;
    }

    double angle, radius = 1.0/3.0;
    double dtheta = 5.0*M_PI/180.0;
 
    for( int i = 0; i < 72; i++) {
        outfile << 5*i << "  ";
        ostringstream oss;
        oss << "anim";
        if( i < 10) oss << "0";
        oss << i << ".vtk";
        JMeshIO::saveAs( mesh, oss.str());
        patchTest(mesh);

        for( int j = 0; j < 4; j++) {
            xyz = nodes[4+j]->getXYZCoords();
            angle = atan2(xyz[1] - 0.5, xyz[0]-0.5) + dtheta;
            xyz[0] = 0.5 + radius*cos(angle);
            xyz[1] = 0.5 + radius*sin(angle);
            nodes[4+j]->setXYZCoords( xyz );
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void quadPatch2()
{
    JMeshPtr mesh = JMesh::newObject();

    JNodeSequence nodes(8);
    for( int i = 0; i < 8; i++) {
        nodes[i] = JVertex::newObject();
        nodes[i]->setID(i);
    }
    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    nodes[0]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    nodes[1]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;
    nodes[2]->setXYZCoords(xyz);

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;
    nodes[3]->setXYZCoords(xyz);

    for( int i = 0; i < 4; i++) {
        double t = 1.7*M_PI + i*0.5*M_PI;
        xyz[0] = 0.5 + 0.25*cos(t);
        xyz[1] = 0.5 + 0.25*sin(t);
        xyz[2] = 0.0;
        nodes[i+4]->setXYZCoords(xyz);
    }

    mesh->addObjects(nodes);

    JFaceSequence faces(5);
    faces[0] = Quadrilateral::newObject( nodes[0], nodes[1], nodes[5], nodes[4] );
    faces[1] = Quadrilateral::newObject( nodes[1], nodes[2], nodes[6], nodes[5] );
    faces[2] = Quadrilateral::newObject( nodes[2], nodes[3], nodes[7], nodes[6] );
    faces[3] = Quadrilateral::newObject( nodes[3], nodes[0], nodes[4], nodes[7] );
    faces[4] = Quadrilateral::newObject( nodes[4], nodes[5], nodes[6], nodes[7] );
    mesh->addObjects( faces);

    JElasticity2D  elastic;

    elastic.setOrder(elemOrder);
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle(0);
    elastic.checkConvexity(0);
    elastic.setMesh(mesh);

    JNodeSequence boundnodes;
    mesh->getTopology()->getBoundary(boundnodes);

    assert( !boundnodes.empty() ) ;

    double u, v;
    for( int i = 0; i < boundnodes.size(); i++) {
        const Point3D &xyz =  boundnodes[i]->getXYZCoords();
        affineFunc(xyz, u,v);
        elastic.fixX( boundnodes[i], u);
        elastic.fixY( boundnodes[i], v);
    }
    elastic.solve();

    vector<double> usol, vsol;
    elastic.getDisplacements(usol, vsol);

    size_t numnodes = mesh->getSize(0);
    cout << fixed << setprecision(15);
    vector<double> error( numnodes);
    double maxerror = 0.0;
    for( int i = 0; i < numnodes; i++) {
        const Point3D &xyz =  mesh->getNodeAt(i)->getXYZCoords();
        affineFunc(xyz, u,v);
        double du = fabs(usol[i] - u)/( fabs(u) + 1.0E-15);
        double dv = fabs(vsol[i] - v)/( fabs(v) + 1.0E-15);
        error[i] = sqrt( du*du + dv*dv);
        maxerror = max(maxerror, du);
        maxerror = max(maxerror, dv);
    }
    cout << "Maximum Error : " << 100.0*maxerror << endl;

    cout << "Info: Done and saving mesh in file : model1.vtk" << endl;
    elastic.saveAs("model1.vtk", "V");
    mesh->deleteAll();
}
/////////////////////////////////////////////////////////////////////////////////////
void quadPatch3()
{
    JMeshPtr mesh = JMesh::newObject();

    JNodeSequence nodes(9);
    for( int i = 0; i < 9; i++) {
        nodes[i] = JVertex::newObject();
        nodes[i]->setID(i);
    }
    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    nodes[0]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    nodes[1]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;
    nodes[2]->setXYZCoords(xyz);

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;
    nodes[3]->setXYZCoords(xyz);

    xyz[0] = 0.5;
    xyz[1] = 0.5;
    xyz[2] = 0.0;
    nodes[4]->setXYZCoords(xyz);

    xyz[0] = 0.20;
    xyz[1] = 0.5;
    xyz[2] = 0.0;
    nodes[5]->setXYZCoords(xyz);

    xyz[0] = 0.50;
    xyz[1] = 0.20;
    xyz[2] = 0.0;
    nodes[6]->setXYZCoords(xyz);

    xyz[0] = 0.80;
    xyz[1] = 0.50;
    xyz[2] = 0.0;
    nodes[7]->setXYZCoords(xyz);

    xyz[0] = 0.50;
    xyz[1] = 0.80;
    xyz[2] = 0.0;
    nodes[8]->setXYZCoords(xyz);

    mesh->addObjects(nodes);

    JFaceSequence faces(6);
    faces[0] = Quadrilateral::newObject( nodes[0], nodes[1], nodes[4], nodes[6] );
    faces[1] = Quadrilateral::newObject( nodes[1], nodes[2], nodes[7], nodes[4] );
    faces[2] = Quadrilateral::newObject( nodes[2], nodes[3], nodes[4], nodes[8] );
    faces[3] = Quadrilateral::newObject( nodes[3], nodes[0], nodes[5], nodes[4] );
    faces[4] = Quadrilateral::newObject( nodes[0], nodes[6], nodes[4], nodes[5] );
    faces[5] = Quadrilateral::newObject( nodes[2], nodes[8], nodes[4], nodes[7] );
    mesh->addObjects( faces);

    JElasticity2D  elastic;

    elastic.setOrder(elemOrder);
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle(0);
    elastic.checkConvexity(0);
    elastic.setMesh(mesh);

    JNodeSequence boundnodes;
    mesh->getTopology()->getBoundary(boundnodes);

    assert( !boundnodes.empty() ) ;

    double u, v;
    for( int i = 0; i < boundnodes.size(); i++) {
        const Point3D &xyz =  boundnodes[i]->getXYZCoords();
        affineFunc(xyz, u,v);
        elastic.fixX( boundnodes[i], u);
        elastic.fixY( boundnodes[i], v);
    }
    elastic.solve();

    vector<double> usol, vsol;
    elastic.getDisplacements(usol, vsol);

    size_t numnodes = mesh->getSize(0);
    cout << fixed << setprecision(15);
    vector<double> error( numnodes);
    double maxerror = 0.0;
    for( int i = 0; i < numnodes; i++) {
        const Point3D &xyz =  mesh->getNodeAt(i)->getXYZCoords();
        affineFunc(xyz, u,v);
        double du = fabs(usol[i] - u)/( fabs(u) + 1.0E-15);
        double dv = fabs(vsol[i] - v)/( fabs(v) + 1.0E-15);
        error[i] = sqrt( du*du + dv*dv);
        maxerror = max( maxerror, du);
        maxerror = max( maxerror, dv);
    }
    cout << "Maximum error : " << 100*maxerror << endl;

    cout << "Info: Done and saving mesh in file : model1.vtk" << endl;
    elastic.saveAs("model1.vtk", "V");
    mesh->deleteAll();

}
/////////////////////////////////////////////////////////////////////////////////////
void patchTest()
{
    int  xlength = 1.0;
    int  ylength = 1.0;

    ofstream ofile( "test.poly", ios::out);

    ofile << "4  2  0  0" << endl;
    ofile << " 0  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << " 1  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << " 2  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << " 3  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << "4  1 " << endl;
    ofile << " 0  0  1  1 "  << endl;
    ofile << " 1  1  2  2 "  << endl;
    ofile << " 2  2  3  3 "  << endl;
    ofile << " 3  3  0  4 "  << endl;
    ofile << " 0 " << endl;
    ofile.close();

    system( "triangle -peq30a0.1 test.poly");

    JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");

    if( meshType == 4) {
        JMeshPtr quadmesh= AllQuadMeshGenerator::SimpleTris2Quads(mesh);
        mesh = quadmesh;
    }

    int etype = mesh->getTopology()->isHomogeneous(2);
    if( etype == 3)
        cout << "Info: Solving problem with a triangle mesh " << endl;

    if( etype == 4)
        cout << "Info: Solving problem with a quad mesh " << endl;

    JElasticity2D  elastic;

    elastic.setOrder(elemOrder);
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle(0);
    elastic.setMesh(mesh);

    JNodeSequence boundnodes;
    mesh->getTopology()->getBoundary(boundnodes);

    assert( !boundnodes.empty() ) ;

    double u, v;
    for( int i = 0; i < boundnodes.size(); i++) {
        const Point3D &xyz =  boundnodes[i]->getXYZCoords();
        affineFunc(xyz, u,v);
        elastic.fixX( boundnodes[i], u);
        elastic.fixY( boundnodes[i], v);
    }
    elastic.solve();

    vector<double> usol, vsol;
    elastic.getDisplacements(usol, vsol);

    size_t numnodes = mesh->getSize(0);
    cout << fixed << setprecision(15);
    vector<double> error( numnodes);
    for( int i = 0; i < numnodes; i++) {
        const Point3D &xyz =  mesh->getNodeAt(i)->getXYZCoords();
        affineFunc(xyz, u,v);
        double du = usol[i] - u;
        double dv = vsol[i] - v;
        error[i] = sqrt( du*du + dv*dv);
        cout << xyz[0] << "  " << xyz[1] << "  " << error[i] << endl;
    }
    double maxerror = *max_element( error.begin(), error.end() );
    cout << "Maximum error in patch test" << maxerror << endl;

    cout << "Info: Done and saving mesh in file : model1.vtk" << endl;
    elastic.saveAs("model1.vtk", "V");
    mesh->deleteAll();
}

////////////////////////////////////////////////////////////////////////////////

int simpleQuadTangle()
{
    JMeshPtr mesh;
    JMeshPtr trimesh = JMeshIO::readFile( "./FEMTest/tangle1.xml");

    JNodePtr v6 = trimesh->getNodeAt(6);
    JNodePtr v7 = trimesh->getNodeAt(7);

    Point3D xyz;

    mesh = AllQuadMeshGenerator::SimpleTris2Quads(trimesh);
    mesh->getTopology()->search_boundary();
    if( randomTangle) {
        Point3D &p6 = v6->getXYZCoords();
        Point3D &p7 = v7->getXYZCoords();
        for( int i = 0; i < 50; i++) {
            p6[0] = 10.0 + 1.77*i;
            p7[0] = 90.0 - 1.77*i;
            outfile << p6[0]/100.0 << " ";
            ostringstream oss;
            oss << "anim";
            if( i < 10) oss << "0";
            oss << i << ".vtk";
            JMeshIO::saveAs(mesh, oss.str());
            patchTest(mesh);
        }
    } else {
        patchTest(mesh);
        JMeshIO::saveAs(mesh, "model.vtk");
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////

int simpleTriTangle()
{
    JMeshPtr mesh = JMeshIO::readFile( "./FEMTest/tangle1.xml");

    JNodePtr v6 = mesh->getNodeAt(6);
    JNodePtr v7 = mesh->getNodeAt(7);

    cout << "#Primary Gauss Points (1,3,7)" << endl;
    cin >> numGauss1;

    cout << "#Secondary Gauss Points (1,3,7)" << endl;
    cin >> numGauss2;

    Point3D xyz;
    xyz[2]  = 0.0;
    for( int k = 0; k < 95; k++) {
        xyz[0] = 1.0+ 1.03*k;
        xyz[1] = 50.0;
        outfile  << xyz[0]/100.0 << "  ";
        v6->setXYZCoords(xyz);
        xyz[0] = 99- 1.03*k;
        xyz[1] = 50.0;
        v7->setXYZCoords(xyz);
        patchTest(mesh);
        ostringstream oss;
        oss << "anim";
        if( k < 10)  oss << "0";
        if( k < 100) oss << "0";
        oss << k << ".vtk";
        JMeshIO::saveAs(mesh, oss.str());
    }
    mesh->deleteAll();
    exit(0);

    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

int testRandomTangle()
{
    JNodeSequence boundnodes;

    int    grid_dim[] = {20, 20};
    double length[]   = {1, 1};
    double origin[]   = {0.0, 0.0};

    cout << "Give grid size (mxn)" << endl;
    cin >> grid_dim[0] >> grid_dim[1];

    JMeshPtr mesh= AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);

    if( meshType == 3)  {
        JMeshPtr trimesh = AllTriMeshGenerator::fromQuadMesh(mesh,4);
        mesh = trimesh;
    }

    int etype = mesh->getTopology()->isHomogeneous(2);
    if( etype == 3) {
        cout << "Info: Solving problem with a triangle mesh " << endl;
        cout << "#Primary  Gauss Points (1,3,7)" << endl;
        cin >> numGauss1;
    }

    if( etype == 4) {
        cout << "Info: Solving problem with a quad  mesh " << endl;
        cout << "#Primary  Gauss Points (4,9)" << endl;
        cin >> numGauss1;
    }

    if( checkTangle ) {
        cout << "#Secondary Gauss Points (1,3,7)" << endl;
        cin >> numGauss2;
    }

    int numnodes = mesh->getSize(0);

    if( randomTangle) {
        double xc = 0.5*length[0];
        double yc = 0.5*length[1];
        double rcut = length[0]/5.0;
        int nCount = 0;
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0] - xc;
            double dy = xyz[1] - yc;
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = xc + (-rcut + 2.0*rcut*drand48());
                xyz[1] = yc + (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
                nCount++;
            }
//         if( nCount == 10) break;
        }
    }
    JMeshIO::saveAs(mesh, "model.vtk");
    patchTest(mesh);
    exit(0);
}

//////////////////////////////////////////////////////////////////////////////////
void BeamVal( double x, double y, double &u, double &v)
{

}
//////////////////////////////////////////////////////////////////////////////////

int testBendingBeam()
{
    int grid_dim[]   = {5, 4};    //  Mediaum
    double length[]  = {1, 0.2};
    double origin[]  = {0.0, 0.0};
    origin[0] = -0.5*length[0];
    origin[1] = -0.5*length[1];

    int gridresol;
    grid_dim[0] = 40;
    grid_dim[1] = 10;
/*
    cout << "Grid resolution (1: Coarse 2: medium 3: fine)" << endl;
    cin  >> gridresol;
    switch(gridresol) {
    case 1:
        grid_dim[0] = 11;
        grid_dim[1] = 7;
        break;
    case 2:
        grid_dim[0] = 20;
        grid_dim[1] = 8;
        break;
    case 3:
        grid_dim[0] = 100;
        grid_dim[1] = 40;
        break;
    }
*/

    JMeshPtr quadmesh= AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);
    JMeshPtr mesh = quadmesh;

    if( meshType == 3) {
        JMeshPtr trimesh = AllTriMeshGenerator::fromQuadMesh(quadmesh,4);
        mesh = trimesh;
    }

/*
    int fid;
    cout << "Which face you want to split" << endl;
    cin >> fid;
    if( fid >= 0) SplitFace( mesh, mesh->getFaceAt(fid));
*/

    int numnodes = mesh->getSize(0);

    // Xshift = -40.0
    // yshift =  10.0

    int loc;
    double xc = 0.0;
    double yc = 0.0;
    if( randomTangle) {
        cout << "Location of tangle elements" << endl;
        cout << "0  :  center      " << endl;
        cout << "1  :  lower left  " << endl;
        cout << "2  :  lower right " << endl;
        cout << "3  :  upper right " << endl;
        cout << "4  :  upper left  " << endl;
        cin >> loc;
        switch(loc) {
        case 1:
            xc = -40.0;
            yc = -10.0;
            break;
        case 2:
            xc =  40.0;
            yc = -10.0;
            break;
        case 3:
            xc =  40.0;
            yc =  10.0;
            break;
        case 4:
            xc = -40.0;
            yc =  10.0;
            break;
        }
        double rcut = 8.0;
        int nCount = 0;
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0] - xc;
            double dy = xyz[1] - yc;
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = xc + (-rcut + 2.0*rcut*drand48());
                xyz[1] = yc + (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
                nCount++;
            }
//            if( nCount == 10) break;
        }
    }
    cout << "tangled mesh stored in model.vtk" << endl;

    JMeshIO::saveAs(mesh, "model.vtk");

    JEdgeSequence edges;

    JElasticity2D  elastic;
    elastic.setOrder(1);
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle(checkTangle);
//  elastic.checkConvexity(0);
    elastic.setMesh(mesh);
    elastic.usePositiveJacobian(absJacobian);

    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
//    elastic.setYForce(edges, 100.0);
     elastic.setXForce(edges, 1.0);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.fixEdges(edges);

    elastic.setShapeFamily(shapeFamily);

    elastic.solve();
    elastic.computeStress();

    vector<double> usol, vsol;
    elastic.getDisplacements(usol, vsol);

    cout << setprecision(10) << endl;

    int id;
    double sval;

    id = grid_dim[0]*grid_dim[1] - grid_dim[0];
    elastic.getValueAt( 'S', id, sval);
    cout << id << "  Stress Value at point D:  " <<  sval << endl;

    id = grid_dim[0]*ceil(0.5*grid_dim[1]) - 1;
    double uval, vval;
    elastic.getValueAt( 'U', id, uval);
    elastic.getValueAt( 'V', id, vval);
    cout << id << " Displacment values at point C:  " <<  uval << "  " << vval << endl;

    cout << "New Data stored in Udata.vtk" << endl;
    elastic.saveAs("Udata.vtk", "U");

    cout << "New Data stored in Vdata.vtk" << endl;
    elastic.saveAs("Vdata.vtk", "V");

    mesh->deleteAll();
    exit(0);

#ifdef LATER
    JLocallyInjectiveMap localmap;
    localmap.setMesh(mesh);
    localmap.setBoundaryPreservation(0);

    JNodeSequence nodes;

    Array3D trans;
    trans[0] = 0.0;
    trans[1] = 0.0;
    trans[2] = 0.0;

    val  = 1;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 3;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    val  = 5;
    int index = 0;
    double sign = 1.0;

    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 90; i++) {
            mopt.improveShapes();
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);
            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }

    sign = -1.0;
    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 90; i++) {
            mopt.improveShapes();
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);
            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }
#endif

}
//////////////////////////////////////////////////////////////////////////////

int testCircle()
{
    double  radius   = 10.0;

    int  npoints = 150;

    ofstream ofile( "test.poly", ios::out);

    ofile << npoints  << " 2  0  0" << endl;

    double dtheta =  2.0*M_PI/(double)npoints;

    for( int i = 0; i < npoints; i++)  {
        double x = radius*cos(i*dtheta);
        double y = radius*sin(i*dtheta);
        ofile << i  << " " << x << "  " << y << endl;
    }

    ofile << npoints <<  "  1 " << endl;
    for( int i = 0; i < npoints; i++)
        ofile << i << "  " << i << " " << (i+1)%npoints << " 1 " << endl;

    ofile << " 0 " << endl;
    ofile.close();
    int err = system( "triangle -peq30a2.1 test.poly");
    assert(!err);

    JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");

    if( randomTangle) {
        double rcut = 4.0;
        int numnodes = mesh->getSize(0);
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0];
            double dy = xyz[1];
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = (-rcut + 2.0*rcut*drand48());
                xyz[1] = (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
            }
        }
    }
    JMeshIO::saveAs(mesh, "tmp1.vtk");
    JMeshNonlinearOptimization mopt;
    mopt.setBoundaryPreservation(1);
    mopt.setMesh(mesh);
    mopt.untangle();
    JMeshIO::saveAs(mesh, "tmp2.vtk");
}

/////////////////////////////////////////////////////////////////////////////////////

int testSquareCircle1()
{
    JElasticity2D  elastic;
    elastic.setOrder(1);

    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);

    JMeshNonlinearOptimization mopt;
    mopt.setBoundaryPreservation(1);

    double  xlength  = 100.0;
    double  ylength  = 40.0;
    double  radius   = 10.0;

    int  npoints = 100;
    double dtheta = 2*M_PI/( double)npoints;

    ofstream ofile( "test.poly", ios::out);

    ofile << npoints+4 << " 2  0  0" << endl;

    for( int i = 0; i < npoints; i++)  {
        double x = radius*cos(i*dtheta);
        double y = radius*sin(i*dtheta);
        ofile << i  << " " << x << "  " << y << endl;
    }

    ofile << npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << npoints + 4 <<  "  1 " << endl;
    for( int i = 0; i < npoints; i++)
        ofile << i << "  " << i << " " << (i+1)%npoints << " 5 " << endl;

    ofile << npoints   << "  " << npoints   << "  " << npoints+1 <<  "  1 " << endl;
    ofile << npoints+1 << "  " << npoints+1 << "  " << npoints+2 <<  "  2 " << endl;
    ofile << npoints+2 << "  " << npoints+2 << "  " << npoints+3 <<  "  3 " << endl;
    ofile << npoints+3 << "  " << npoints+3 << "  " << npoints   <<  "  4 " << endl;

    int landmark2 = npoints + 2;

    ofile << " 1 " << endl;
    ofile << " 0  0.0 0.0 " << endl;
    ofile.close();
    system( "triangle -peq30a5.0 test.poly");

    JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");

    int numnodes = mesh->getSize(0);
    assert( numnodes);

    double xc, yc;
    if( randomTangle) {
        cout << "Random tangle position " << endl;
        cout << "1  :   top of the cylinder " << endl;
        cout << "2  :   bottom of the cylinder " << endl;
        cout << "3  :   left and right of the cylinder " << endl;
        int pos;
        cin >> pos;
        switch(pos) {
        case 1:
            xc = 0.0;
            yc = 15.0;
            break;
        case 2:
            xc = 0.0;
            yc =-15.0;
            break;
        case 3:
            xc = -30.0;
            yc =   0.0;
            break;
        }

        double rcut = 3.0;
        int nCount = 0;
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0] - xc;
            double dy = xyz[1] - yc;
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = xc + (-rcut + 2.0*rcut*drand48());
                xyz[1] = yc + (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
                nCount++;
            }
        }

        if( pos == 3 ) {
            xc =  30.0;
            for( int i = 0; i < numnodes; i++) {
                JNodePtr v = mesh->getNodeAt(i);
                Point3D xyz = v->getXYZCoords();
                double dx = xyz[0] - xc;
                double dy = xyz[1] - yc;
                double r = sqrt(dx*dx + dy*dy);
                if( r < rcut)  {
                    xyz[0] = xc + (-rcut + 2.0*rcut*drand48());
                    xyz[1] = yc + (-rcut + 2.0*rcut*drand48());
                    v->setXYZCoords(xyz);
                    nCount++;
                }
            }
        }
    }

    if( meshType == 4) {
        JMeshPtr quadmesh= AllQuadMeshGenerator::SimpleTris2Quads(mesh);
        mesh = quadmesh;
    }

    int landmark1 = 0, shift = 0;
    double startAngle = shift*dtheta;
    Point3D xyz;
    double mindist = std::numeric_limits<double>::max();
    for( int i = 0; i < npoints; i++)  {
        xyz[0] = radius*cos(startAngle + i*dtheta);
        xyz[1] = radius*sin(startAngle + i*dtheta);
        xyz[2] = 0.0;
        double dx = xyz[0];
        double dy = xyz[1] - radius;
        double dl = sqrt(dx*dx + dy*dy);
        if( dl < mindist) {
            landmark1 = i;
            mindist = dl;
        }
        mesh->getNodeAt(i)->setXYZCoords(xyz);
    }

    elastic.checkTangle(checkTangle);
    elastic.setMesh(mesh);
    elastic.usePositiveJacobian(absJacobian);

    JEdgeSequence edges;

    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.setXForce(edges, 100.0);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.fixEdges(edges);

    elastic.solve();
    elastic.computeStress();

    cout << "Storing data in model.vtk" << endl;
    elastic.saveAs("model.vtk", "S");

    cout << setprecision(10);

    xyz = mesh->getNodeAt(landmark1)->getXYZCoords();
    cout << xyz[0] << "  " << xyz[1] << endl;

    double sval;
    elastic.getValueAt( 'S', landmark1, sval);
    cout << "Stress Value at point D:  " <<  sval << endl;

    xyz = mesh->getNodeAt(landmark2)->getXYZCoords();
    cout << xyz[0] << "  " << xyz[1] << endl;

    double uval, vval;
    elastic.getValueAt( 'U', landmark2, uval);
    elastic.getValueAt( 'V', landmark2, vval);
    cout << "Displacment values at point C:  " <<  uval << "  " << vval << endl;

}
//////////////////////////////////////////////////////////////////////////////
int LinearElasticitySquareCircle2( JMeshPtr &mesh)
{
    double exact = 10.1987;
    JElasticity2D  elastic;
    int etype = mesh->getTopology()->isHomogeneous(2);
    elastic.setOrder(elemOrder);
    elastic.setElementType(etype);

    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle( checkTangle );
    elastic.checkConvexity(0);
    elastic.setNumPrimaryFaceGaussPoints(numGauss1);
    elastic.setNumSecondaryFaceGaussPoints(numGauss2);
    elastic.setMesh(mesh);

    int val;
    JEdgeSequence edges;

    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.setXForce(edges, 100.0);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.fixEdges(edges);

    elastic.solve();

    elastic.computeStress();

    cout << "Storing data in model.vtk" << endl;
    elastic.saveAs("model.vtk", "S");

    val  = 5;
    mesh->getEntities("Boundary", val, edges);
    JNodeSequence nodes;
    JMeshTopology::getEntitySet( edges, nodes);

    double ymax = 0.0;
    int    landmark = 0;
    Point3D xyz;
    for( size_t i = 0; i < nodes.size(); i++) {
        xyz = nodes[i]->getXYZCoords();
        if( xyz[1] > ymax ) {
            ymax = xyz[1];
            landmark = nodes[i]->getID();
        }
    }
    cout << "Landmark node " << landmark << endl;
    xyz =  mesh->getNodeAt(landmark)->getXYZCoords();
    cout << "Landmark position " << xyz[0] << " " << xyz[1] << endl;

    double sval;
    elastic.getValueAt( 'S', landmark, sval);
    cout << "Stress Value at point D:  " <<  sval << endl;
    outfile << 100.0*fabs(sval - exact)/exact << endl;
//  outfile << sval << endl;
}
//////////////////////////////////////////////////////////////////////////////

int testSquareCircle2()
{
    JMeshPtr trimesh;

    JEdgeSequence edges;
    JNodeSequence nodes;

    int    grid_dim[] = {4, 2};
    double length[]   = {100, 40};
    double origin[]   = {0.0, 0.0};

    JMeshPtr mesh= AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);

    cout << "#Primary  Gauss Points (4,9)" << endl;
    cin >> numGauss1;

    JEdgePtr edge;
    JFacePtr face;
    int   val;

//  First cell: Three sides exposed to boundary
    face = mesh->getFaceAt(0);

    val = 1;
    edge= face->getEdgeAt(0);
    face->setAttribute("Boundary", val);

    val = 3;
    edge = face->getEdgeAt(2);
    edge->setAttribute("Boundary", val);

    val = 4;
    edge = face->getEdgeAt(3);
    edge->setAttribute("Boundary", val);

//  Second Cell: Two side exposed to boundary

    face = mesh->getFaceAt(1);
    val  = 1;
    edge = face->getEdgeAt(0);
    edge->setAttribute("Boundary", val);

    val  = 3;
    edge = face->getEdgeAt(2);
    edge->setAttribute("Boundary", val);

//  Third cell: Three sides exposed to boundary
    face = mesh->getFaceAt(2);

    val = 1;
    edge= face->getEdgeAt(0);
    face->setAttribute("Boundary", val);

    val = 2;
    edge = face->getEdgeAt(1);
    edge->setAttribute("Boundary", val);

    val = 3;
    edge = face->getEdgeAt(2);
    edge->setAttribute("Boundary", val);

    QuadRefiner refiner;
    refiner.setMesh(mesh);
    JNodeSequence newnodes;
    JFaceSequence newfaces;

    edge = mesh->getFaceAt(1)->getEdgeAt(0);

    int bmark;
    int err = edge->getAttribute("Boundary", bmark);
    assert( err == 0);

    refiner.refine15( mesh->getFaceAt(1), newnodes, newfaces);

    edge = newfaces[0]->getEdgeAt(0);
    err = edge->getAttribute("Boundary", bmark);
    assert( err == 0);

    val = 5;
    for( int i = 0; i < 4; i++) {
        edge = newfaces[4]->getEdgeAt(i);
        edge->setAttribute("Boundary", val);
    }

    assert( newfaces.size() == 5);

    newfaces[4]->setStatus(JMeshEntity::REMOVE);
    mesh->pruneFaces();
    grid_dim[0] = 11;
    grid_dim[1] = 11;
    refiner.refineAll( grid_dim );

    if( meshType == 3)  {
        trimesh = AllTriMeshGenerator::fromQuadMesh(mesh,4);
        mesh = trimesh;
    }

    JLocallyInjectiveMap localmap;
    if( meshType == 4)  {
        trimesh = AllTriMeshGenerator::fromQuadMesh(mesh,4);
        localmap.setMesh(trimesh);
    } else
        localmap.setMesh(mesh);

    localmap.setBoundaryPreservation(0);

    Array3D trans;
    trans[0] = 0.0;
    trans[1] = 0.0;
    trans[2] = 0.0;
    cout << "Exit now : CSV " << endl;
    exit(0);

    /*
        for( int j = 0; j < 4; j++) {
            mesh->getEntities("Boundary", j+1, edges);
            JMeshTopology::getEntitySet( edges, nodes);
            for( int i = 0; i < nodes.size(); i++) {
                localmap.translate(nodes[i], trans);
            }
        }
    */

    val  = 5;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);

    double xc = 0.0;
    double yc = 0.0;
    int nsize = nodes.size();
    for( int i = 0; i < nsize; i++) {
        JNodePtr v = nodes[i];
        Point3D xyz = v->getXYZCoords();
        xc += xyz[0];
        yc += xyz[1];
    }

    xc = xc/(double)nsize;
    yc = yc/(double)nsize;
    double radius = 10.0;

    circleCenter[0] = xc;
    circleCenter[1] = yc;
    circleRadius    = radius;

//   double start_angle = 2.0*M_PI/(double)nodes.size();
    double start_angle = rotAngle*M_PI/180.0;

    int    pos = 0;
    double mindist = std::numeric_limits<double>::max();
    double t, dx, dy, dl;
    Point3D xyz;

    for( int i = 0; i < nodes.size(); i++) {
        JNodePtr v = nodes[i];
        xyz = v->getXYZCoords();
        dx = xyz[0] - xc;
        dy = xyz[1] - yc;

        t = atan2(dy,dx);
        xyz[0] = xc + radius*cos(t + start_angle);
        xyz[1] = yc + radius*sin(t + start_angle);
        xyz[2] = 0.0;
        v->setAttribute("TargetPos", xyz);

        dx = xyz[0] - xc ;
        dy = xyz[1] - yc - radius;
        dl = sqrt(dx*dx + dy*dy);
        if( dl < mindist) {
            pos = i;
            mindist = dl;
        }
    }
    int landmark = nodes[pos]->getID();
    xyz =  nodes[pos]->getXYZCoords();

    localmap.setMaxIterations(10000);
    localmap.solve();

    JMeshNonlinearOptimization mopt;
    mopt.setBoundaryPreservation(1);
    mopt.setNumIterations(100);
    mopt.setMesh(mesh);
    mopt.improveShapes();

    for( size_t i = 0; i < nodes.size(); i++) {
        JNodePtr v = nodes[i];
        Point3D xyz = v->getXYZCoords();
        double dx = xyz[0] - xc;
        double dy = xyz[1] - yc;
        double r  = sqrt(dx*dx + dy*dy);
        assert( fabs(r - radius) < 1.0E-06);
    }

    circleRadius = radius;

#ifdef CSV
    if( randomTangle) {
        double rcut = 10.0;
        int nCount = 0;
        int numnodes = mesh->getSize(0);
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double xt = 15.0;
            double yt = 20.0;
            double dx = xyz[0] - xt;
            double dy = xyz[1] - yt;
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = xt + (-rcut + 2.0*rcut*drand48());
                xyz[1] = yt + (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
                nCount++;
            }
//         if( nCount == 10) break;
        }
    }


    if( randomTangle ) {
        val  = 5;
        mesh->getEntities("Boundary", val, edges);

        QuadChord qChord;
        qChord.setMesh(mesh);
        qChord.setSeed( edges[0] );
        JFaceSequence qFaces = qChord.getMeshFaces();

        qChord.setSeed( qFaces[3]->getEdgeAt(1) );
        qFaces = qChord.getMeshFaces();

        JEdgeSet eset;
        for( int i = 0; i < qFaces.size(); i++) {
            eset.insert( qFaces[i]->getEdgeAt(1));
            eset.insert( qFaces[i]->getEdgeAt(3));
        }
        map<JNodePtr, Point3D> orgPoints;

        foreach_( JEdgePtr edge, eset) {
            JNodePtr v0 = edge->getNodeAt(0);
            JNodePtr v1 = edge->getNodeAt(1);
            orgPoints[v0] = v0->getXYZCoords();
            orgPoints[v1] = v1->getXYZCoords();
        }

        cout << "#Secondary Gauss Points (1,3,7)" << endl;
        cin >> numGauss2;

        double dr = 0.043;
        Point3D p0, p1;

        for( int i = 0; i < 25; i++) {
            double r = i*dr;
            outfile << r << " ";
            foreach_(JEdgePtr edge, eset) {
                JNodePtr v0 = edge->getNodeAt(0);
                JNodePtr v1 = edge->getNodeAt(1);

                p0[0]   = (1.0-r)*orgPoints[v0][0] + r*orgPoints[v1][0];
                p0[1]   = (1.0-r)*orgPoints[v0][1] + r*orgPoints[v1][1];
                p0[2]   = 0.0;

                p1[0]   = (1.0-r)*orgPoints[v1][0] + r*orgPoints[v0][0];
                p1[1]   = (1.0-r)*orgPoints[v1][1] + r*orgPoints[v0][1];
                p1[2]   = 0.0;

                v0->setXYZCoords(p0);
                v1->setXYZCoords(p1);
            }
//          patchTest(mesh);
            LinearElasticitySquareCircle2( mesh);
            ostringstream oss;
            oss << "anim";
            if( i < 10 ) oss << "0";
            oss << i << ".vtk";
            mesh->saveAs(oss.str());
        }
        exit(0);
    }
    LinearElasticitySquareCircle2( mesh);
//  patchTest(mesh);
#endif
    exit(0);
}
//////////////////////////////////////////////////////////////////////////////


int testSquareCircles3()
{
    double  xlength  = 100.0;
    double  ylength  = 40.0;
    double  radius1  = 10;
    double  radius2  = 5;

    int  npoints = 100;
    double dtheta = 2*M_PI/( double)npoints;

    ofstream ofile( "test.poly", ios::out);

    ofile << 3*npoints+4 << " 2  0  0" << endl;

    int index = 0;
    for( int i = 0; i < npoints; i++)  {
        double x = radius1*cos(i*dtheta);
        double y = radius1*sin(i*dtheta);
        ofile << index++  << " " << x << "  " << y << endl;
    }

    for( int i = 0; i < npoints; i++)  {
        double x = -20 + radius2*cos(i*dtheta);
        double y = radius2*sin(i*dtheta);
        ofile << index++  << " " << x << "  " << y << endl;
    }

    for( int i = 0; i < npoints; i++)  {
        double x = 20 + radius2*cos(i*dtheta);
        double y = radius2*sin(i*dtheta);
        ofile << index++  << " " << x << "  " << y << endl;
    }

    ofile << 3*npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << 3*npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << 3*npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << 3*npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << 3*npoints + 4 <<  "  1 " << endl;
    int offset = 0;
    for( int i = 0; i < npoints; i++)
        ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 5 " << endl;

    offset = npoints;
    for( int i = 0; i < npoints; i++)
        ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 6 " << endl;

    offset = 2*npoints;
    for( int i = 0; i < npoints; i++)
        ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 7 " << endl;

    ofile << 3*npoints   << "  " << 3*npoints   << "  " << 3*npoints+1 <<  "  1 " << endl;
    ofile << 3*npoints+1 << "  " << 3*npoints+1 << "  " << 3*npoints+2 <<  "  2 " << endl;
    ofile << 3*npoints+2 << "  " << 3*npoints+2 << "  " << 3*npoints+3 <<  "  3 " << endl;
    ofile << 3*npoints+3 << "  " << 3*npoints+3 << "  " << 3*npoints   <<  "  4 " << endl;

    ofile << " 3 " << endl;
    ofile << " 0   0.0  0.0 " << endl;
    ofile << " 0  -20.0 0.0 " << endl;
    ofile << " 0   20 0.0.0 " << endl;

    ofile.close();
    system( "triangle -peq30a0.2 test.poly");

    JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");

    int numnodes = mesh->getSize(0);

    if( randomTangle) {
        double rotangle = 10.0*(M_PI/180.0);
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0] + 20;
            double dy = xyz[1];
            double angle = atan2(dy,dx) + rotangle;
            double r = sqrt(dx*dx + dy*dy);
            if( r > 5.1 && r < 7.5 )  {
                double rnew = 5.5 + drand48();
                xyz[0] = rnew*cos(angle) - 20.0;
                xyz[1] = rnew*sin(angle);
                v->setXYZCoords(xyz);
            }
            if( r > 5.5 && r < 6.5 )  {
                double rnew = 5.1 + drand48();
                xyz[0] = rnew*cos(angle) - 20.0;
                xyz[1] = rnew*sin(angle);
                v->setXYZCoords(xyz);
            }
        }

        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0] - 20;
            double dy = xyz[1];
            double angle = atan2(dy,dx) + rotangle;
            double r = sqrt(dx*dx + dy*dy);
            if( r > 5.1 && r < 7.5 )  {
                double rnew = 5.5 + drand48();
                xyz[0] = rnew*cos(angle) + 20.0;
                xyz[1] = rnew*sin(angle);
                v->setXYZCoords(xyz);
            }
            if( r > 5.5 && r < 6.0 )  {
                double rnew = 5.1 + drand48();
                xyz[0] = rnew*cos(angle) + 20.0;
                xyz[1] = rnew*sin(angle);
                v->setXYZCoords(xyz);
            }
        }
    }

    if( meshType == 4) {
        JMeshPtr quadmesh= AllQuadMeshGenerator::SimpleTris2Quads(mesh);
        mesh = quadmesh;
    }

    JElasticity2D  elastic;
    elastic.setOrder(1);
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
    elastic.checkTangle(checkTangle);
    elastic.setMesh(mesh);

    JNodeSequence nodes;
    JEdgeSequence edges;

    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.setXForce(edges, 100.0);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.fixEdges(edges);

    elastic.solve();
    elastic.computeStress();
    cout << "Data stored in model.vtk" << endl;
    elastic.saveAs("model.vtk", "S");

    exit(0);

    JLocallyInjectiveMap localmap;
    localmap.setMesh(mesh);
    localmap.setBoundaryPreservation(0);


    Array3D trans;
    trans[0] = 0.0;
    trans[1] = 0.0;
    trans[2] = 0.0;

    val  = 1;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 3;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 5;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    val  = 5;
    double sign = 1.0;

    index = 0;
    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 90; i++) {
            mopt.improveShapes();

            val  = 6;
            mesh->getEntities("Boundary", val, edges);
            JMeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);

            val  = 7;
            mesh->getEntities("Boundary", val, edges);
            JMeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);

            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }

    sign = -1.0;
    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 90; i++) {
            mopt.improveShapes();

            val  = 6;
            mesh->getEntities("Boundary", val, edges);
            JMeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);

            val  = 7;
            mesh->getEntities("Boundary", val, edges);
            JMeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);

            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }


    /*
        JElasticity2D  elastic;
        elastic.setOrder(1);
        elastic.setYoungModulus(100);
        elastic.setPoissonRatio(0.25);
        elastic.checkTangle(0);
        elastic.setMesh(mesh);

        JEdgeSequence edges;

        int val;
        val  = 2;
        mesh->getEntities("Boundary", val, edges);
        elastic.setXForce(edges, 100.0);

        cout << edges.size() << endl;

        val  = 4;
        mesh->getEntities("Boundary", val, edges);
        elastic.setXForce(edges, -100.0);

        val  = 5;
        mesh->getEntities("Boundary", val, edges);
        elastic.fixEdges(edges);

        elastic.solve();
        elastic.saveAs("anim.vtk");
        return mesh;
    */
}
//////////////////////////////////////////////////////////////////////////////

#ifdef CSV

int  npoints = 100;
double dtheta = 2*M_PI/( double)npoints;

ofstream ofile( "test.poly", ios::out);

ofile << 3*npoints+4 << " 2  0  0" << endl;

int index = 0;
for( int i = 0; i < npoints; i++)
{
    double x = radius1*cos(i*dtheta);
    double y = radius1*sin(i*dtheta);
    ofile << index++  << " " << x << "  " << y << endl;
}

for( int i = 0; i < npoints; i++)
{
    double x = -20 + radius2*cos(i*dtheta);
    double y = radius2*sin(i*dtheta);
    ofile << index++  << " " << x << "  " << y << endl;
}

for( int i = 0; i < npoints; i++)
{
    double x = 20 + radius2*cos(i*dtheta);
    double y = radius2*sin(i*dtheta);
    ofile << index++  << " " << x << "  " << y << endl;
}

ofile << 3*npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
ofile << 3*npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
ofile << 3*npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
ofile << 3*npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

ofile << 3*npoints + 4 <<  "  1 " << endl;
int offset = 0;
for( int i = 0; i < npoints; i++)
    ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 5 " << endl;

offset = npoints;
for( int i = 0; i < npoints; i++)
    ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 6 " << endl;

offset = 2*npoints;
for( int i = 0; i < npoints; i++)
    ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 7 " << endl;

ofile << 3*npoints   << "  " << 3*npoints   << "  " << 3*npoints+1 <<  "  1 " << endl;
ofile << 3*npoints+1 << "  " << 3*npoints+1 << "  " << 3*npoints+2 <<  "  2 " << endl;
ofile << 3*npoints+2 << "  " << 3*npoints+2 << "  " << 3*npoints+3 <<  "  3 " << endl;
ofile << 3*npoints+3 << "  " << 3*npoints+3 << "  " << 3*npoints   <<  "  4 " << endl;

ofile << " 3 " << endl;
ofile << " 0   0.0  0.0 " << endl;
ofile << " 0  -20.0 0.0 " << endl;
ofile << " 0   20 0.0.0 " << endl;

ofile.close();
system( "triangle -peq30a0.2 test.poly");

Mesh *mesh = Mesh::newObject();
mesh->readFromFile( "test.1.ele");

int numnodes = mesh->getSize(0);

if( randomTangle)
{
    double rotangle = 10.0*(M_PI/180.0);
    for( int i = 0; i < numnodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        Point3D xyz = v->getXYZCoords();
        double dx = xyz[0] + 20;
        double dy = xyz[1];
        double angle = atan2(dy,dx) + rotangle;
        double r = sqrt(dx*dx + dy*dy);
        if( r > 5.1 && r < 7.5 )  {
            double rnew = 5.5 + drand48();
            xyz[0] = rnew*cos(angle) - 20.0;
            xyz[1] = rnew*sin(angle);
            v->setXYZCoords(xyz);
        }
        if( r > 5.5 && r < 6.5 )  {
            double rnew = 5.1 + drand48();
            xyz[0] = rnew*cos(angle) - 20.0;
            xyz[1] = rnew*sin(angle);
            v->setXYZCoords(xyz);
        }
    }

    for( int i = 0; i < numnodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        Point3D xyz = v->getXYZCoords();
        double dx = xyz[0] - 20;
        double dy = xyz[1];
        double angle = atan2(dy,dx) + rotangle;
        double r = sqrt(dx*dx + dy*dy);
        if( r > 5.1 && r < 7.5 )  {
            double rnew = 5.5 + drand48();
            xyz[0] = rnew*cos(angle) + 20.0;
            xyz[1] = rnew*sin(angle);
            v->setXYZCoords(xyz);
        }
        if( r > 5.5 && r < 6.0 )  {
            double rnew = 5.1 + drand48();
            xyz[0] = rnew*cos(angle) + 20.0;
            xyz[1] = rnew*sin(angle);
            v->setXYZCoords(xyz);
        }
    }
}

if( meshType == 4)
{
    Mesh *quadmesh= AllQuadMeshGenerator::SimpleTris2Quads(mesh);
    mesh = quadmesh;
}

JElasticity2D  elastic;
elastic.setOrder(1);
    elastic.setYoungModulus(youngModulus);
    elastic.setPoissonRatio(poissonRatio);
elastic.checkTangle(checkTangle);
elastic.setMesh(mesh);

JNodeSequence nodes;
JEdgeSequence edges;

int val;
val  = 2;
mesh->getEntities("Boundary", val, edges);
assert( !edges.empty() ) ;
elastic.setXForce(edges, 100.0);

val  = 4;
mesh->getEntities("Boundary", val, edges);
assert( !edges.empty() ) ;
elastic.fixEdges(edges);

elastic.solve();
elastic.computeStress();
cout << "Data stored in model.vtk" << endl;
elastic.saveAs("model.vtk", "S");

exit(0);

JLocallyInjectiveMap localmap;
localmap.setMesh(mesh);
localmap.setBoundaryPreservation(0);


Array3D trans;
trans[0] = 0.0;
trans[1] = 0.0;
trans[2] = 0.0;

val  = 1;
mesh->getEntities("Boundary", val, edges);
MeshTopology::getEntitySet( edges, nodes);
localmap.translate(nodes, trans);

val  = 2;
mesh->getEntities("Boundary", val, edges);
MeshTopology::getEntitySet( edges, nodes);
localmap.translate(nodes, trans);

val  = 3;
mesh->getEntities("Boundary", val, edges);
MeshTopology::getEntitySet( edges, nodes);
localmap.translate(nodes, trans);

val  = 4;
mesh->getEntities("Boundary", val, edges);
MeshTopology::getEntitySet( edges, nodes);
localmap.translate(nodes, trans);

val  = 5;
mesh->getEntities("Boundary", val, edges);
MeshTopology::getEntitySet( edges, nodes);
localmap.translate(nodes, trans);

JMeshNonlinearOptimization mopt;
mopt.setMesh(mesh);
mopt.setBoundaryPreservation(1);

val  = 5;
double sign = 1.0;

index = 0;
for( int j = 0; j < 2; j++)
{
    for( int i = 0; i < 90; i++) {
        mopt.improveShapes();

        val  = 6;
        mesh->getEntities("Boundary", val, edges);
        MeshTopology::getEntitySet( edges, nodes);
        localmap.rotate(nodes, 2.0*sign);

        val  = 7;
        mesh->getEntities("Boundary", val, edges);
        MeshTopology::getEntitySet( edges, nodes);
        localmap.rotate(nodes, 2.0*sign);

        double maxDist = localmap.solve();
        cout << "MaxDist  " << maxDist << endl;
        ostringstream oss;
        oss << "./animateData/anim";
        if( index < 10)  oss << "0";
        if( index < 100) oss << "0";
        oss << index << ".vtk";
        mesh->saveAs( oss.str());
        index++;
    }
    sign *= -1.0;
}

sign = -1.0;
for( int j = 0; j < 2; j++)
{
    for( int i = 0; i < 90; i++) {
        mopt.improveShapes();

        val  = 6;
        mesh->getEntities("Boundary", val, edges);
        MeshTopology::getEntitySet( edges, nodes);
        localmap.rotate(nodes, 2.0*sign);

        val  = 7;
        mesh->getEntities("Boundary", val, edges);
        MeshTopology::getEntitySet( edges, nodes);
        localmap.rotate(nodes, 2.0*sign);

        double maxDist = localmap.solve();
        cout << "MaxDist  " << maxDist << endl;
        ostringstream oss;
        oss << "./animateData/anim";
        if( index < 10)  oss << "0";
        if( index < 100) oss << "0";
        oss << index << ".vtk";
        mesh->saveAs( oss.str());
        index++;
    }
    sign *= -1.0;
}


/*
    JElasticity2D  elastic;
    elastic.setOrder(1);
    elastic.setYoungModulus(100);
    elastic.setPoissonRatio(0.25);
    elastic.checkTangle(0);
    elastic.setMesh(mesh);

    JEdgeSequence edges;

    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    elastic.setXForce(edges, 100.0);

    cout << edges.size() << endl;

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    elastic.setXForce(edges, -100.0);

    val  = 5;
    mesh->getEntities("Boundary", val, edges);
    elastic.fixEdges(edges);

    elastic.solve();
    elastic.saveAs("anim.vtk");
    return mesh;
*/
}
//////////////////////////////////////////////////////////////////////////////


int testmodel7()
{
    int  xlength = 100.0;
    int  ylength = 50.0;
    double  radius  = 9.5;

    int  npoints = 150;
    double dtheta = 2*M_PI/( double)npoints;

    ofstream ofile( "test.poly", ios::out);

    ofile << 2*npoints+4 << " 2  0  0" << endl;

    int index = 0;
    for( int i = 0; i < npoints; i++)  {
        double x = -30.0 + radius*cos(i*dtheta);
        double y =  10.0 + radius*sin(i*dtheta);
        ofile << index++  << " " << x << "  " << y << endl;
    }

    for( int i = 0; i < npoints; i++)  {
        double x =  30.0 + radius*cos(i*dtheta);
        double y = -10.0 + radius*sin(i*dtheta);
        ofile << index++  << " " << x << "  " << y << endl;
    }

    ofile << 2*npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << 2*npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << 2*npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << 2*npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << 2*npoints + 4 <<  "  1 " << endl;
    int offset = 0;
    for( int i = 0; i < npoints; i++)
        ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 5 " << endl;

    offset = npoints;
    for( int i = 0; i < npoints; i++)
        ofile << offset+i << "  " << offset+i << " " << offset+(i+1)%npoints << " 6 " << endl;

    ofile << 2*npoints   << "  " << 2*npoints   << "  " << 2*npoints+1 <<  "  1 " << endl;
    ofile << 2*npoints+1 << "  " << 2*npoints+1 << "  " << 2*npoints+2 <<  "  2 " << endl;
    ofile << 2*npoints+2 << "  " << 2*npoints+2 << "  " << 2*npoints+3 <<  "  3 " << endl;
    ofile << 2*npoints+3 << "  " << 2*npoints+3 << "  " << 2*npoints   <<  "  4 " << endl;

    ofile << " 2 " << endl;
    ofile << " 0   -30.0  10.0 " << endl;
    ofile << " 0    30.0 -10.0 " << endl;

    ofile.close();
    system( "triangle -peq30a1.0 test.poly");

    Mesh *mesh = Mesh::newObject();
    mesh->readFromFile( "test.1.ele");

    JLocallyInjectiveMap localmap;
    localmap.setMesh(mesh);
    localmap.setBoundaryPreservation(0);

    JNodeSequence nodes;
    JEdgeSequence edges;

    Array3D trans;
    trans[0] = 0.0;
    trans[1] = 0.0;
    trans[2] = 0.0;

    int val  = 1;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 3;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    val  = 5;
    double sign = 1.0;

    index = 0;
    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 30; i++) {
            mopt.improveShapes();

            val = 5;
            trans[0] = 1.0;
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.translate(nodes, trans);

            val = 6;
            trans[0] = -1.0;
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.translate(nodes, trans);

            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }

    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 30; i++) {
            mopt.improveShapes();

            val = 5;
            trans[0] = -1.0;
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.translate(nodes, trans);

            val = 6;
            trans[0] =  1.0;
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.translate(nodes, trans);

            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }


    /*
        sign = -1.0;
        for( int j = 0; j < 2; j++) {
            for( int i = 0; i < 90; i++) {
                mopt.improveShapes();
                mesh->getEntities("Boundary", val, edges);
                MeshTopology::getEntitySet( edges, nodes);
                localmap.rotate(nodes, 2.0*sign);
                double maxDist = localmap.solve();
                cout << "MaxDist  " << maxDist << endl;
                ostringstream oss;
                oss << "./animateData/anim";
                if( index < 10)  oss << "0";
                if( index < 100) oss << "0";
                oss << index << ".vtk";
                mesh->saveAs( oss.str());
                index++;
            }
            sign *= -1.0;
        }
    */

    /*
        JElasticity2D  elastic;
        elastic.setOrder(1);
        elastic.setYoungModulus(100);
        elastic.setPoissonRatio(0.25);
        elastic.checkTangle(0);
        elastic.setMesh(mesh);

        JEdgeSequence edges;

        int val;
        val  = 2;
        mesh->getEntities("Boundary", val, edges);
        elastic.setXForce(edges, 100.0);

        cout << edges.size() << endl;

        val  = 4;
        mesh->getEntities("Boundary", val, edges);
        elastic.setXForce(edges, -100.0);

        val  = 5;
        mesh->getEntities("Boundary", val, edges);
        elastic.fixEdges(edges);

        elastic.solve();
        elastic.saveAs("anim.vtk");
        return mesh;
    i*/
}
//////////////////////////////////////////////////////////////////////////////

void testmodel8()
{
    int  xlength = 40.0;
    int  ylength = 100.0;
    double  radius  = 10.0;
    double  height = 30.0;

    int  npoints = 100;
    double dtheta = 2*M_PI/( double)npoints;

    ofstream ofile( "test.poly", ios::out);

    ofile << npoints+4 << " 2  0  0" << endl;

    for( int i = 0; i < npoints; i++)  {
        double x = radius*cos(i*dtheta);
        double y = height+radius*sin(i*dtheta);
        ofile << i  << " " << x << "  " << y << endl;
    }
    ofile << npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << npoints + 4 <<  "  1 " << endl;
    for( int i = 0; i < npoints; i++)
        ofile << i << "  " << i << " " << (i+1)%npoints << " 5 " << endl;

    ofile << npoints   << "  " << npoints   << "  " << npoints+1 <<  "  1 " << endl;
    ofile << npoints+1 << "  " << npoints+1 << "  " << npoints+2 <<  "  2 " << endl;
    ofile << npoints+2 << "  " << npoints+2 << "  " << npoints+3 <<  "  3 " << endl;
    ofile << npoints+3 << "  " << npoints+3 << "  " << npoints   <<  "  4 " << endl;

    ofile << " 1 " << endl;
    ofile << " 0  0.0 " << height << endl;
    ofile.close();
    system( "triangle -peq30a10.0 test.poly");

    Mesh *mesh = Mesh::newObject();
    mesh->readFromFile( "test.1.ele");
    mesh->saveAs("anim.vtk");
    mesh->saveAs("FEMTest/model8.xml");
    exit(0);

    JElasticity2D  elastic;
    elastic.setOrder(1);
    elastic.setYoungModulus(100);
    elastic.setPoissonRatio(0.25);
    elastic.checkTangle(0);
    elastic.setMesh(mesh);

    JEdgeSequence edges;

    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    elastic.setXForce(edges, 100.0);

    cout << edges.size() << endl;

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    elastic.setXForce(edges, -100.0);

    val  = 5;
    mesh->getEntities("Boundary", val, edges);
    elastic.fixEdges(edges);

    elastic.solve();
    elastic.saveAs("anim.vtk", "U");
}

////////////////////////////////////////////////////////////////////////////////
void testmodel9()
{
    vector<size_t> permute;

    int grid_dim[]   = {20, 10};
    double length[]  = {10, 5};
    double origin[]  = {0.0, 0.0};
    /*
        Mesh * mesh = getStructuredMesh(2, grid_dim, length, origin);
        mesh->saveAs( "tmp.vtk");
        exit(0);
    */
    Mesh *mesh = Mesh::newObject();
    mesh->readFromFile( "./FEMTest/model9.xml");

    mesh->getTopology()->search_boundary();

    int numnodes = mesh->getSize(0);
    permute.resize(numnodes);
    for( int i = 0; i < numnodes; i++)
        permute[i] = i;
    SwapPositions( mesh->getNodeAt(116), mesh->getNodeAt(95), permute);
    mesh->saveAs( "tmp.vtk");

    JElasticity2D  elastic;
    elastic.setOrder(1);
    elastic.setYoungModulus(100);
    elastic.setPoissonRatio(0.25);
    elastic.setMesh(mesh);
    elastic.checkTangle(0);

    JNodeSequence boundnodes;
    mesh->getTopology()->getBoundary( boundnodes);

    for( int i = 0; i < boundnodes.size(); i++) {
        const Point3D &xyz =  boundnodes[i]->getXYZCoords();
        double x = xyz[0];
        double y = xyz[1];
        double u =  0.1*x - 0.2*y;
        double v = -0.3*x + 0.5*y;
        elastic.fixX( boundnodes[i], u);
        elastic.fixY( boundnodes[i], v);
    }
    elastic.solve();

    vector<double> usol, vsol;
    elastic.getDisplacements(usol, vsol);
    double umax = 0.0;
    double vmax = 0.0;
    double errmax = 0.0;


    numnodes = mesh->getSize(0);
    cout << fixed << setprecision(15);
    vector<double> error( numnodes);
    for( int i = 0; i < numnodes; i++) {
        const Point3D &xyz =  mesh->getNodeAt(i)->getXYZCoords();
        double x = xyz[0];
        double y = xyz[1];
        double u =  0.1*x - 0.2*y;
        double v = -0.3*x + 0.5*y;
        double du = usol[i] - u;
        double dv = vsol[i] - v;
        error[i] = sqrt( du*du + dv*dv);
        /*
        //      cout << i << "  " << usol[i] << "  " << vsol[i] << " Exact " << u << "  " << v << endl;
                umax = max( umax, fabs(usol[i] - u));
                vmax = max( vmax, fabs(vsol[i] - v));
        */
    }

    if( umax > 1.0E-10  || vmax > 1.0E-10)
        cout << "Value do not match" <<  endl;

    cout << "Error " << endl;
    sort( error.begin(), error.end() ) ;
    for( int i = 0; i < error.size(); i++)
        cout << i << "  " << error[i] << endl;


    elastic.saveAs("tmp.vtk", "U");
}

/////////////////////////////////////////////////////////////////////////////////////////////

int testmodel10()
{
    Mesh *mesh = Mesh::newObject();

    vector<size_t> permute;

    /*
        int grid_dim[]   = {41, 21};
        double length[]  = {10, 5};
        double origin[]  = {0.0, 0.0};
        mesh = getStructuredMesh(2, grid_dim, length, origin);
        mesh->saveAs( "tmp.vtk");
        exit(0);
    */

    mesh->readFromFile( "./FEMTest/model10.xml");

    JEdgeSequence edges;

    permute.resize(mesh->getSize(0));

    JElasticity2D  elastic;
    elastic.setOrder(1);
    elastic.setYoungModulus(100);
    elastic.setPoissonRatio(0.25);
    elastic.checkTangle(0);
    elastic.setMesh(mesh);

    int landmarks[] = { 0, 10, 20, 30, 40,
                        245, 450, 491, 860,
                        850, 840, 830, 820,
                        451, 410, 205
                      };


    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    elastic.setYForce(edges, 100.0);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    elastic.fixEdges(edges);

    mesh->getTopology()->search_boundary();

    Point3D xyz;
    for( int i = 0; i < 100; i++) {
        JNodePtr vtx = mesh->getRandomNode();
        if( !vtx->isBoundary() )  {
            xyz[0] = + 8*drand48();
            xyz[1] = + 4*drand48();
            xyz[2] =  0.0;
            vtx->setXYZCoords(xyz);
        }
    }

    elastic.solve();
    elastic.saveAs("anim.vtk", "S");
    cout << "New Data stored in anim.vtk" << endl;

    cout << "Values at landmarks" << endl;
    double nodeval;
    for( int i = 0; i < 16; i++) {
        xyz = mesh->getNodeAt(landmarks[i])->getXYZCoords();
        cout << xyz[0] << "  " << xyz[1] <<  "  ";
        elastic.getValueAt( 'U', landmarks[i], nodeval);
        cout << nodeval << "  ";
        elastic.getValueAt( 'V', landmarks[i], nodeval);
        cout << nodeval << "  ";
        elastic.getValueAt( 'S', landmarks[i], nodeval);
        cout << nodeval << "  ";
        cout << endl;
    }

    exit(0);

    JLocallyInjectiveMap localmap;
    localmap.setMesh(mesh);
    localmap.setBoundaryPreservation(0);

    JNodeSequence nodes;

    Array3D trans;
    trans[0] = 0.0;
    trans[1] = 0.0;
    trans[2] = 0.0;

    val  = 1;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 3;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    MeshTopology::getEntitySet( edges, nodes);
    localmap.translate(nodes, trans);

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    val  = 5;
    int index = 0;
    double sign = 1.0;

    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 90; i++) {
            mopt.improveShapes();
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);
            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }

    sign = -1.0;
    for( int j = 0; j < 2; j++) {
        for( int i = 0; i < 90; i++) {
            mopt.improveShapes();
            mesh->getEntities("Boundary", val, edges);
            MeshTopology::getEntitySet( edges, nodes);
            localmap.rotate(nodes, 2.0*sign);
            double maxDist = localmap.solve();
            cout << "MaxDist  " << maxDist << endl;
            ostringstream oss;
            oss << "./animateData/anim";
            if( index < 10)  oss << "0";
            if( index < 100) oss << "0";
            oss << index << ".vtk";
            mesh->saveAs( oss.str());
            index++;
        }
        sign *= -1.0;
    }
}
#endif

///////////////////////////////////////////////////////////////////////////////////////
int testSimpleTangle()
{
    if( meshType == 3) simpleTriTangle();
    if( meshType == 4) simpleQuadTangle();
}

///////////////////////////////////////////////////////////////////////////////////////
int testSelfTanglePatch()
{
    
/*
   int    grid_dim[] = {5, 4};
   double length[]   = {1, 1};
   double origin[]   = {0.0, 0.0};

   cout << "Give grid size (mxn)" << endl;
   cin >> grid_dim[0] >> grid_dim[1];
   JMeshPtr mesh= AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);
*/
 
    JMeshPtr mesh = JMesh::newObject();
    JNodeSequence nodes(5);
    nodes[0] = JVertex::newObject();
    nodes[1] = JVertex::newObject();
    nodes[2] = JVertex::newObject();
    nodes[3] = JVertex::newObject();
    nodes[4] = JVertex::newObject();

    Point3D p3d;
    p3d[2] = 0.0;

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    nodes[0]->setXYZCoords(p3d);

    p3d[0] = 1.0;
    p3d[1] = 0.0;
    nodes[1]->setXYZCoords(p3d);

    p3d[0] = 1.0;
    p3d[1] = 1.0;
    nodes[2]->setXYZCoords(p3d);

    p3d[0] = 0.0;
    p3d[1] = 1.0;
    nodes[3]->setXYZCoords(p3d);

    p3d[0] = 0.1;
    p3d[1] = 0.9;

    mesh->addObjects(nodes);
    mesh->addObject( Quadrilateral::newObject(nodes[0], nodes[1], nodes[4], nodes[3] ));
    mesh->addObject( Quadrilateral::newObject(nodes[1], nodes[2], nodes[3], nodes[4] ));

    int nx = 5;
    int ny = 5;
    double dx = 1.0/(double)nx;
    double dy = 1.0/(double)ny;

    for( int j = 1; j < ny-1; j++)  {
        for( int i = 1; i < nx-1; i++)  {
             p3d[0]  = i*dx;
             p3d[1]  = j*dy;
             nodes[4]->setXYZCoords(p3d);
             patchTest(mesh);
         }
    }
    
    JMeshIO::saveAs(mesh, "model.vtk");
    return 0;
}


//////////////////////////////////////////////////////////////////////////////////
void inputParams()
{
    outfile.open("error.dat", ios::out);

    ifstream infile( "test.config", ios::in);

    if( infile.fail() ) {
    cout << "Select Shape familty : (0) All classical (1) Mixed ( classical + barycentri) (2) All barycentric " << endl;
    cin >> shapeFamily;

    while(1) {
    cout << "MeshType : (Tri=3, Quad=4 ) " << endl;
    cin >> meshType;
    if( meshType == 3 || meshType == 4 ) break;
    }

    while(1) {
    cout << "Shape order : (Linear = 1: Quadratric 2" << endl;
    cin >> elemOrder;
    if( elemOrder == 1 || elemOrder == 2) break;
    }

    cout << "Create random tangle ( Yes = 1, No = 0)" << endl;
    cin  >> randomTangle;

    cout << "Check mesh tangle : (Yes=1 , No = 0)" << endl;
    cin >> checkTangle;

    cout << "Do you want to use traditional positive Jacobian(yes = 1, no =0)" << endl;
    cin >> absJacobian;

    cout << "Field Type)" << endl;
    cin >> fieldType;

    return;
    }

    int tgauss, qgauss;

    string str;
    while( !infile.eof() ) {
          infile >> str;
          if( str == "FieldType"    ) infile >>  fieldType;
          if( str == "MeshType"     ) infile >>  meshType;
          if( str == "TriGauss1"    ) infile >>  tgauss;
          if( str == "TriGauss2"    ) infile >>  numGauss2;
          if( str == "QuadGauss"    ) infile >>  qgauss;
          if( str == "AbsJacobian"  ) infile >>  absJacobian;
          if( str == "ShapeFunc"    ) infile >>  shapeFamily;
          if( str == "ElemOrder"   )  infile >>  elemOrder;
          if( str == "RandomTangle" ) infile >>  randomTangle;
          if( str == "CheckTangle"  ) infile >>  checkTangle;
     }

     if( meshType == 3) 
         numGauss1 = tgauss;
     else
         numGauss1 = qgauss;
     
}
//////////////////////////////////////////////////////////////////////////////////
#endif

int main()
{
/*
    JFacePtr f = Quadrilateral::getCanonical();

    exit(0);

    cout << "TestModel " << endl;
    cout << " 1 : Self Tangle test   " << endl;
    cout << " 2 : Simple tangle       " << endl;
    cout << " 3 : Random tangle       " << endl;
    cout << " 4 : Rotate Quad         " << endl;
    cout << " 5 : Linear elasticity (Bending Beam )  " << endl;
    cout << " 6 : Circle               " << endl;
    cout << " 7 : Linear elasticity (square+ 1 circle)  " << endl;
    cout << " 8 : Linear elasticity (square+ 2 circle)  " << endl;
    cin >> testmodel;

    inputParams();
    

    switch( testmodel ) {
    case 1:
        testSelfTanglePatch();
        break;
    case 2:
        testSimpleTangle();
        break;
    case 3:
        testRandomTangle();
        break;
    case 4:
        testRotateQuadPatch1();
        break;
    case 5:
        testBendingBeam();
        break;
    case 6:
        testCircle();
        break;
    case 7:
        testSquareCircle1();
        break;
    case 8:
        testSquareCircles3();
        break;
    }
*/

    return 0;
}
