#include "Mesh.hpp"
#include "DelaunayMesh.hpp"
#include "MeshLaplacian.hpp"

using namespace Jaal;

int main( int argc, char **argv)
{
     DelaunayMesh delmesh;
     vector<Point2D>  points;

     int npoints = 1000;

     points.resize(npoints);
     for( int i = 0; i < npoints; i++) {
          points[i][0] = drand48();
          points[i][1] = drand48();
     }
     Mesh *mesh = delmesh.addPoints(points);
     mesh->saveAs("tmp.xml");

     JLaplaceMeshSmoother msmooth;
     msmooth.setMesh(mesh);
     msmooth.setNumIterations(1000);
     msmooth.setBoundaryPreservation(1);
     msmooth.smoothAll();
     mesh->saveAs("tmp1.xml");
     exit(0);

     JMeshNonlinearOptimization mopt;
     mopt.setMesh(mesh);
     mopt.setNumIterations(1000);
     mopt.improveShapes();
     mesh->saveAs("tmp1.xml");

     vector<double>  coords1;
     vector<size_t>  l2g;

     mesh->getGeometry()->getCoordsArray(coords1, l2g);
     delmesh.retriangulate(mesh);
     mesh->saveAs("tmp2.xml");
}
