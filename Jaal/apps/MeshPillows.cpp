#include "Poisson2D.hpp"
#include <iomanip>
#include"MeshOptimization.hpp"
#include "MeshRefine.hpp"
#include "MeshLaplacian.hpp"
#include "BinaryTreeMatch.hpp"

/////////////////////////////////////////////////////////////////////////////////

int testCircle()
{
    double  radius = 1.0;

    int  npoints = 8;
    double dtheta = 2*M_PI/( double)npoints;

    JNodeSequence nodes(npoints+1);

    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    Vertex *v  = Vertex::newObject();
    v->setID(0);
    v->setXYZCoords(xyz);
    nodes[0] = v;
    int bmark;
    for( int i = 0; i < npoints; i++)  {
        xyz[0] = radius*cos(i*dtheta);
        xyz[1] = radius*sin(i*dtheta);
        Vertex *v  = Vertex::newObject();
        v->setXYZCoords(xyz);
        v->setID(i+1);
        nodes[i+1] = v;
        v->setAttribute("Constraint", bmark);
    }
    Mesh *mesh = Mesh::newObject();
    mesh->addObjects( nodes);

    JFaceSequence faces(npoints);
    for( int i = 0; i < npoints; i++) {
        Vertex *v0 = nodes[0];
        Vertex *v1 = nodes[(i+1)];
        Vertex *v2 = nodes[(i+1)%npoints+1];
        faces[i]   = Triangle::newObject( v0, v1, v2);
    }
    mesh->addObjects(faces);

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    JMeshLaplaceSmoother laplace;
    laplace.setMesh(mesh);
    laplace.setNumIterations(2);

    QuadRefiner refiner(mesh);
    for( int i = 0; i < 50; i++) {
        refiner.insert_boundary_pillows();
        laplace.smooth();
        mopt.improveShapes();
    }

    mesh->saveAs("tmp.vtk");
}
////////////////////////////////////////////////////////////////////////////////

int testSqrCircle()
{
    vector<size_t> permute;

    double  xlength  = 50.0;
    double  ylength  = 50.0;
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
        ofile << i << "  " << i << " " << (i+1)%npoints << " 1 " << endl;

    ofile << npoints   << "  " << npoints   << "  " << npoints+1 <<  " 2 " << endl;
    ofile << npoints+1 << "  " << npoints+1 << "  " << npoints+2 <<  " 2 " << endl;
    ofile << npoints+2 << "  " << npoints+2 << "  " << npoints+3 <<  " 2 " << endl;
    ofile << npoints+3 << "  " << npoints+3 << "  " << npoints   <<  " 2 " << endl;

    ofile << " 1 " << endl;
    ofile << " 0  0.0 0.0 " << endl;
    ofile.close();
    system( "triangle -peq30a50.0 test.poly");

    Mesh *mesh = Mesh::newObject();
    mesh->readFromFile( "test.1.ele");

    JBinaryTreeMatch tmatch;
    Mesh *quadmesh = tmatch.getQuadMesh(mesh);

    quadmesh->saveAs("tmp.vtk");
    exit(0);

    





    mesh->getTopology()->search_boundary();

    JNodeSequence bnodes;
    mesh->getTopology()->getBoundary(bnodes);
    int bmark = 1;
    for( Vertex *vtx: bnodes) vtx->setAttribute("Constraint", bmark);

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    JMeshLaplaceSmoother laplace;
    laplace.setMesh(mesh);
    laplace.setNumIterations(100);
    QuadRefiner refiner(mesh);
    for( int i = 0; i < 9; i++) {
        refiner.insert_boundary_pillows();
        laplace.smooth();
        mopt.improveShapes();
    }

    mesh->saveAs("tmp.vtk");

}

/////////////////////////////////////////////////////////////////////////////////////

int main()
{
    testSqrCircle();

    return 0;
}

