#include <iomanip>

#include "Poisson2D.hpp"
#include "AllQuadMeshGenerator.hpp"
/*
#include"MeshOptimization.hpp"
#include "MeshRefine.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AllQuadMeshGenerator.hpp"
*/


struct JConstantField : public JFieldFunction
{
    double  getScalar( const Point2D &) const {
        return  1.0;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
int testSquare()
{
    int grid_dim[]   = {8, 8};
    double length[]  = {1, 1};
    double origin[]  = {-1.0, -1.0};
    JMeshPtr qmesh = AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);
//  JMeshPtr  mesh = AllTriMeshGenerator::fromQuadMesh( qmesh, 2);
    JMeshPtr  mesh = qmesh;
    mesh->getTopology()->searchBoundary();
#ifdef CSV

    /*
        cout << "Solving Poisson equation in a square domain " << endl;
        int  xlength = 1.0;
        int  ylength = 1.0;
        int  npoints = 0;

        ofstream ofile( "test.poly", ios::out);

        ofile << "4 2 0 0" << endl;

    //    ofile << npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    //    ofile << npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    //    ofile << npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    //    ofile << npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

        ofile << npoints   << "  " <<  0.0 << " " <<  0.0 << endl;
        ofile << npoints+1 << "  " <<  1.0 << " " <<  0.0 << endl;
        ofile << npoints+2 << "  " <<  1.0 << " " <<  1.0 << endl;
        ofile << npoints+3 << "  " <<  0.0 << " " <<  1.0 << endl;

        ofile << "4 1 " << endl;
        ofile << npoints   << "  " << npoints   << "  " << npoints+1 <<  " 1 " << endl;
        ofile << npoints+1 << "  " << npoints+1 << "  " << npoints+2 <<  " 2 " << endl;
        ofile << npoints+2 << "  " << npoints+2 << "  " << npoints+3 <<  " 3 " << endl;
        ofile << npoints+3 << "  " << npoints+3 << "  " << npoints   <<  " 4 " << endl;

        ofile << " 0 " << endl;
        ofile.close();

        system( "triangle -peq30a0.0001 test.poly");
        JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");
    */

    JPoisson2D  poisson;
    poisson.setOrder(2);
    poisson.setMesh(mesh);
    poisson.checkTangle(0);

    size_t numnodes = mesh->getSize(0);
    for( size_t i  = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isBoundary() ) {
            const Point3D &xyz = vtx->getXYZCoords();
            double x = xyz[0];
            double y = xyz[1];
            poisson.setDirichletValue(vtx,  0.0);
        }
    }

/*
    JEdgeSequence edges;
    double val = 0.0;
    mesh->getEntities("Boundary", 1, edges);
    poisson.setNeumannValue( edges, val);

    val = -4.0;
    mesh->getEntities("Boundary", 3, edges);
    poisson.setNeumannValue( edges, val);
*/
    boost::shared_ptr<JFieldFunction>  field( new JConstantField);
    poisson.setFieldFunction( field );

    poisson.solve();
    poisson.saveAs("poisson.vtk", "U");
    exit(0);

    double ucal, uexact;
    double maxerror = 0.0;
    for( size_t i  = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D  &xyz = vtx->getXYZCoords();
        double x = xyz[0];
        double y = xyz[1];
        vtx->getAttribute("U", ucal);
        uexact =  1 + x*x + 2*y*y;
        maxerror = max( maxerror, fabs(ucal-uexact));
        cout << maxerror << endl;
    }
    cout << "Maximum Error " << maxerror << endl;
#endif

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int testLShape()
{
/*
    cout << "Solving Poisson equation in a square domain " << endl;
    int  xlength = 1.0;
    int  ylength = 1.0;
    int  npoints = 0;

    ofstream ofile( "test.poly", ios::out);

    ofile << "6 2 0 0" << endl;

    ofile << " 0  0.0 1.0 " << endl;
    ofile << " 1  1.0 1.0 " << endl;
    ofile << " 2  1.0 0.0 " << endl;
    ofile << " 3  2.0 0.0 " << endl;
    ofile << " 4  2.0 2.0 " << endl;
    ofile << " 5  0.0 2.0 " << endl;

    ofile << "6 1 " << endl;
    ofile << "1  0 1   1"  << endl;
    ofile << "2  1 2   2" << endl;
    ofile << "3  2 3   3" << endl;
    ofile << "4  3 4   4" << endl;
    ofile << "5  4 5   5" << endl;
    ofile << "6  5 0   6" << endl;
    ofile << " 0 " << endl;
    ofile.close();

    system( "triangle -peq30a0.001 test.poly");
    JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");

    JPoisson2D  poisson;
    poisson.setOrder(1);
    poisson.setMesh(mesh);
    poisson.checkTangle(0);

    JEdgeSequence edges;

    double val = 0.0;
    mesh->getEntities("Boundary", 6, edges);
    assert( !edges.empty() );
    poisson.setDirichletValue( edges, val);

    val = 1.0;
    mesh->getEntities("Boundary", 3, edges);
    assert( !edges.empty() );
    poisson.setDirichletValue( edges, val);

    boost::shared_ptr<JFieldFunction>  field( new JConstantField);
    poisson.setFieldFunction( field );

    poisson.solve();
    poisson.saveAs("poisson.vtk", "U");
*/

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int testFork()
{
#ifdef CSV
    vector<size_t> permute;

    double len1 = 0.5;
    double len2 = 1.5;
    double len3 = 0.5;
    double len4 = 0.5;
    double len5 = 2.0;
    double len6 = 0.5;

    int  npoints = 0;

    ofstream ofile( "test.poly", ios::out);

    ofile << "12 2 0 0" << endl;

    ofile << "0  "   << "  " <<  0.5*len1 << " " <<  "0.0" << endl;
    ofile << "1 "    << "  " <<  0.5*len1 << " " <<  len2  << endl;
    ofile << "2 "   << "  " <<   0.5*len1+len3 << " " <<  len2  << endl;
    ofile << "3 "   << "  " <<   0.5*len1+len3 << " " <<  len2+len4+len5  << endl;
    ofile << "4 "   << "  " <<   0.5*len1+len3-len6  << " " <<  len2+len4+len5  << endl;
    ofile << "5 "   << "  " << 0.5*len1+len3-len6 << " " <<  len2+len4  << endl;

    ofile << "6 "   << "  " << -0.5*len1+len3-len6 << " " <<  len2+len4  << endl;
    ofile << "7 "   << "  " << -(0.5*len1+len3-len6) << " " <<  len2+len4+len5  << endl;
    ofile << "8 "   << "  " << -(0.5*len1+len3) << " " <<  len2+len4+len5  << endl;
    ofile << "9 "   << "  " << -(0.5*len1+len3) << " " <<  len2  << endl;
    ofile << "10 "   << "  " << -0.5*len1 << " " <<  len2  << endl;
    ofile << "11 "   << "  " << -0.5*len1 << " " <<  "0.0" << endl;

    ofile << "12 1 " << endl;
    ofile << "0   11  0  1 " << endl;
    ofile << "1   0   1  2 " << endl;
    ofile << "2   1   2  3 " << endl;
    ofile << "3   2   3  4 " << endl;
    ofile << "4   3   4  5 " << endl;
    ofile << "5   4   5  6 " << endl;
    ofile << "6   5   6  7 " << endl;
    ofile << "7   6   7  8 " << endl;
    ofile << "8   7   8  9 " << endl;
    ofile << "9   8   9 10 " << endl;
    ofile << "10  9  10 11 " << endl;
    ofile << "11  10 11 12 " << endl;

    ofile << " 0 " << endl;
    ofile.close();

    system( "triangle -peq30a0.005 test.poly");

    JMeshPtr mesh = MeshIO::readFile("test.1.ele");


    /*
        JPoisson2D  poisson;
        poisson.setOrder(1);
        poisson.checkTangle(0);
        poisson.setMesh(mesh);

        int val;

        JEdgeSequence edges;
        mesh->getEntities("Boundary", 1, edges);
        assert( !edges.empty() ) ;
        poisson.setField(edges, 0.0);

        mesh->getEntities("Boundary", 5, edges);
        assert( !edges.empty() ) ;
        poisson.setField(edges, 1.0);

        mesh->getEntities("Boundary", 9, edges);
        assert( !edges.empty() ) ;
        poisson.setField(edges, 1.0);

        poisson.solve();

        poisson.solve();
        poisson.saveAs("anim.vtk", "U");
        cout << "New Data stored in anim.vtk" << endl;
    */

#endif

}

//////////////////////////////////////////////////////////////////////////////////

int testSqrCircle()
{
#ifdef CSV
    vector<size_t> permute;

    int  xlength   = 50.0;
    int  ylength   = 50.0;
    double  radius = 10.0;

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
    system( "triangle -peq30a1.0 test.poly");

    JMeshPtr mesh = MeshIO::readFile("test.1.ele");

    /*
        JPoisson2D  poisson;
        poisson.setOrder(1);
        poisson.checkTangle(0);
        poisson.setMesh(mesh);
        poisson.setNumPrimaryFaceGaussPoints(7);

        int val;

        JEdgeSequence edges;
        val  = 1;
        mesh->getEntities("Boundary", val, edges);
        assert( !edges.empty() ) ;
        poisson.setField(edges, 0.0);

        val  = 2;
        mesh->getEntities("Boundary", val, edges);
        assert( !edges.empty() ) ;
        poisson.setField(edges, 100.0);

        poisson.solve();
        poisson.saveAs("anim.vtk", "U");
        cout << "New Data stored in anim.vtk" << endl;
    */
#endif

}

/////////////////////////////////////////////////////////////////////////////////////
void getOverlap()
{
/*
    int grid_dim[]   = {100, 100};
    double length[]  = {2, 2};
    double origin[]  = {-1.0, -1.0};

    JMeshPtr mesh1 = AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);
    JMeshPtr mesh2 = mesh1->deepCopy();

    int nx = grid_dim[0];
    int ny = grid_dim[0];

   double dx = 2.0/(double)(nx-1);
   double dy = 2.0/(double)(ny-1);

   JNodeSequence nodes(4);
   nodes[0] = JVertex::newObject();
   nodes[1] = JVertex::newObject();
   nodes[2] = JVertex::newObject();
   nodes[3] = JVertex::newObject();

   Point3D xyz;
   xyz[2]  = 0.0;

   xyz[0]  = 0.0;
   xyz[1]  = 0.0;
   nodes[0]->setXYZCoords(xyz);

   xyz[0]  = 1.0;
   xyz[1]  = 0.0;
   nodes[1]->setXYZCoords(xyz);

   xyz[0]  =  1.0;
   xyz[1]  =  1.0;
   nodes[2]->setXYZCoords(xyz);

   xyz[0]  = 0.8;
   xyz[1]  = 0.1;
   nodes[3]->setXYZCoords(xyz);
   JFacePtr face = Quadrilateral::newObject( nodes);

   JFEM2D  fem;

   vector<Point2D> uvCoords;
   vector<Array4I> quads;
   fem.splitSelfTangledQuad( face, uvCoords, quads);

   JPoisson2D poisson;
   MatrixXd Ke;
   vector<double> fe;

   cout << " Hello " << endl;
   poisson.integrate_self_tangle( face, Ke, fe);

   double  J;
   Point2D uv;
   Point2D xy;

   double dJ, val;

   size_t numnodes = mesh1->getSize(0);
   for( size_t i = 0; i < numnodes; i++) {
         const JNodePtr &v1 = mesh1->getNodeAt(i);
         const JNodePtr &v2 = mesh2->getNodeAt(i);
         uv = v1->getXYCoords();
         J = fem.getJacobian(face, uv);
         fem.getXYCoords(face, uv, xy);
         val = J < 0.0 ? -1.0 : 1.0;
         v1->setAttribute("Jacobian", val);
         v2->setAttribute("Jacobian", val);
         fem.getXYCoords(face, uv, xy);
         xyz[0] = xy[0];
         xyz[1] = xy[1];
         xyz[2] = 0.0;
         v2->setXYZCoords(xyz);
   }

    JMeshVTKExporter mexp;
    mexp.addNodeAttribute("Jacobian");
    mexp.writeFile(mesh1, "model1.vtk");
    mexp.writeFile(mesh2, "model2.vtk");
*/

}
/////////////////////////////////////////////////////////////////////////////////////
//
int main()
{

#ifdef CSV
    testSquare();
    exit(0);
    vector<Point2D>   xy(4);
    xy[3][0] = 0.0;
    xy[3][1] = 0.0;
    
    xy[0][0] = 1.0;
    xy[0][1] = 0.0;

    xy[1][0] = 0.45;
    xy[1][1] = 0.45;

    xy[2][0] = 0.0;
    xy[2][1] = 1.0;

    ConcaveQuadElement quad;
    int c;
    Point2D ab, uv;
    c = quad.getConcaveCorner( xy);
    quad.getReflectionPoint(xy, ab);

    vector<double> N;

    double a = ab[0];
    double b = ab[1];
    ofstream ofile( "test.poly", ios::out);

    double  U[6], V[6];
    if( c == 0) {
        U[0] = a; U[1] = a;  U[2] = 1;  U[3] = 1; U[4] =-1; U[5] = -1;
        V[0] = b; V[1] = -1; V[2] = -1; V[3] = 1; V[4] = 1;  V[5] = b;
    }

    if( c == 1) {
        U[0] = -1; U[1] =  a; U[2] = a; U[3] = 1; U[4] = 1;  U[5] = -1;
        V[0] = -1; V[1] = -1; V[2] = b; V[3] = b; V[4] = 1;  V[5] =  1;
    }

    if( c == 2) {
        U[0] = -1; U[1] =  1; U[2] = 1; U[3] = a; U[4] = a;  U[5] = -1;
        V[0] = -1; V[1] = -1; V[2] = b; V[3] = b; V[4] = 1;  V[5] =  1;
    }

    if( c == 3) {
        U[0] = -1; U[1] =  1; U[2] = 1; U[3] = a; U[4] = a;  U[5] = -1;
        V[0] = -1; V[1] = -1; V[2] = 1; V[3] = 1; V[4] = b;  V[5] =  b;
    }

    ofile << "6 2 0 0" << endl;
    ofile << " 0  " << U[0] << " " << V[0] << endl;
    ofile << " 1  " << U[1] << " " << V[1] << endl;
    ofile << " 2  " << U[2] << " " << V[2] << endl;
    ofile << " 3  " << U[3] << " " << V[3] << endl;
    ofile << " 4  " << U[4] << " " << V[4] << endl;
    ofile << " 5  " << U[5] << " " << V[5] << endl;

    ofile << "6 1 " << endl;
    ofile << "1  0 1   1"  << endl;
    ofile << "2  1 2   2" << endl;
    ofile << "3  2 3   3" << endl;
    ofile << "4  3 4   4" << endl;
    ofile << "5  4 5   5" << endl;
    ofile << "6  5 0   6" << endl;
    ofile << " 0 " << endl;

    ofile.close();

    system( "triangle -peq30a0.001 test.poly");
    JMeshPtr mesh = JMeshIO::readFile( "test.1.ele");

    for( int i = 0; i < mesh->getSize(0); i++) {
         const JVertexPtr &vtx = mesh->getNodeAt(i);
         uv = vtx->getXYCoords();
         quad.getShapeFunc1(c, ab, uv, N);
         vtx->setAttribute("N0", N[0] );
         vtx->setAttribute("N1", N[1] );
         vtx->setAttribute("N2", N[2] );
         vtx->setAttribute("N3", N[3] );
         vtx->setAttribute("N4", N[4] );
         vtx->setAttribute("N5", N[5] );
         vtx->setAttribute("N6", N[6] );
         vtx->setAttribute("N7", N[7] );
     }

     JMeshVTKExporter vtk0;
     vtk0.addNodeAttribute("N0");
     vtk0.writeFile(mesh, "model0.vtk");

     JMeshVTKExporter vtk1;
     vtk1.addNodeAttribute("N1");
     vtk1.writeFile(mesh, "model1.vtk");

     JMeshVTKExporter vtk2;
     vtk2.addNodeAttribute("N2");
     vtk2.writeFile(mesh, "model2.vtk");

     JMeshVTKExporter vtk3;
     vtk3.addNodeAttribute("N3");
     vtk3.writeFile(mesh, "model3.vtk");

     JMeshVTKExporter vtk4;
     vtk4.addNodeAttribute("N4");
     vtk4.writeFile(mesh, "model4.vtk");

     JMeshVTKExporter vtk5;
     vtk5.addNodeAttribute("N5");
     vtk5.writeFile(mesh, "model5.vtk");

     JMeshVTKExporter vtk6;
     vtk6.addNodeAttribute("N6");
     vtk6.writeFile(mesh, "model6.vtk");

     JMeshVTKExporter vtk7;
     vtk7.addNodeAttribute("N7");
     vtk7.writeFile(mesh, "model7.vtk");
/*
    getOverlap();
    testLShape();
*/
#endif
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////

