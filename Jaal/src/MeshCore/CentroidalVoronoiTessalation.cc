#ifndef CENTOIDVORO_H
#define CENTOIDVORO_H

#include <stdio.h>
#include <stdlib.h>

#include "Mesh.hpp"
#include "SwapTriEdges.hpp"
#include "circumcenter.hpp"
#include "MeshRefine2D.hpp"
#include "MeshPartitioner.hpp"

#include <vector>

using namespace std;
using namespace Jaal;

class CentroidalVoronoiTessalator {
public:
     void smooth( vector<int> &triConnect, vector<double> &coords, int maxiter = 100);
     void smooth( Mesh *m, int iter = 100 );
private:
     struct LVertex {
        Vertex  *vertex;
        void  getNewCoord();
        double setNewCoord() {
             const Point3D &oldpos = vertex->getXYZCoords();
             double dx = fabs( newpos[0] - oldpos[0] );
             double dy = fabs( newpos[1] - oldpos[1] );
             double dz = fabs( newpos[2] - oldpos[2] );
             vertex->setXYZCoords( newpos );
             return Math::max_value( dx, dy, dz );
        }
        JFaceSequence faces;
        Point3D newpos;
     };
     struct EdgeFaces
     {
        Face *faces[2];
        bool contains( const Face *f ) {
             if( faces[0] == f ) return 1;
             if( faces[1] == f)  {
                 std::swap( faces[0], faces[1] );
                 return 1;
             }
             return 0;
             }
     };

     void lloyd_relaxation();
     void build_vertex_rings();

     int maxiter;
     Jaal::Mesh *mesh;
     MeshGeometricOptimization mopt;
     vector<LVertex>  pnodes;
};

void  CentroidalVoronoiTessalator::LVertex::getNewCoord()
{
     Point3D pj, param;
     Array3D angles;

     int nSize = faces.size();
     assert( nSize > 0);

     vector<Point2D> pc( nSize );
     for( int i = 0; i < nSize; i++) {
          Face *f = faces[i];
          assert(f->getSize(0) == 3);
          Point3D p0 = f->getNodeAt(0)->getXYZCoords();
          Point3D p1 = f->getNodeAt(1)->getXYZCoords();
          Point3D p2 = f->getNodeAt(2)->getXYZCoords();
          Math::getTriAngles( p0, p1, p2, angles);
          double maxang = Math::max_value(angles[0], angles[1], angles[2] );
          if( maxang > 90.0) 
              f->getAvgPos( pj );
           else 
             TriCircumCenter2D( &p0[0], &p1[0], &p2[0], &pj[0], &param[0] );
          pc[i][0] = pj[0];
          pc[i][1] = pj[1];
     }
     Point2D c;
     Math::poly_centroid( pc, c );
     newpos[0] = c[0];
     newpos[1] = c[1];
     newpos[2] = 0.0;

}

void CentroidalVoronoiTessalator :: build_vertex_rings()
{
     size_t nSize = mesh->getSize(0);
     int numneighs;

     pnodes.reserve(nSize);
     LVertex lv;
     vector<EdgeFaces> vedgefaces;
     JFaceSequence efaces;
     JNodeSequence vneighs;
     for( size_t i = 0; i < nSize; i++)  {
          Vertex *vertex = mesh->getNodeAt(i);
          if( vertex->isActive() &&  !vertex->isBoundary()  ) {
               lv.vertex = vertex;
               vertex->getRelations( vneighs );
               numneighs = vneighs.size();
               vedgefaces.resize( numneighs );
               for( int j = 0; j < numneighs; j++) {
                    Mesh::getRelations112(vertex, vneighs[j], efaces);
                    assert( efaces.size() == 2 );
                    vedgefaces[j].faces[0] = efaces[0];
                    vedgefaces[j].faces[1] = efaces[1];
               }
               Face *nextface = vedgefaces[0].faces[1];
               for( int j = 0; j  < numneighs-1; j++) {
                    for( int k = j+1; k < numneighs; k++) {
                         if( vedgefaces[k].contains(nextface) ){
                             nextface = vedgefaces[k].faces[1];
                             std::swap( vedgefaces[j], vedgefaces[k] );
                             break;
                         }
                    }
               }
               for( int j = 0; j  < numneighs; j++) {
                   if( vedgefaces[i].faces[1] != vedgefaces[(i+1)%numneighs].faces[0] ) {
                       cout << "Fatal Error: loop not formed: " << endl;
                       exit(0);
                   }
              }

              lv.faces.resize( numneighs );
              for( int j = 0; j  < numneighs; j++) 
                   lv.faces[j] = vedgefaces[j].faces[0];
              pnodes.push_back( lv );
          }
     }
}

void CentroidalVoronoiTessalator :: lloyd_relaxation()
{
     build_vertex_rings();

     double error, maxerror = 0.0;

     size_t nSize = pnodes.size();
     for( int iter = 0; iter < maxiter; iter++) {
     for( size_t i = 0; i < nSize; i++)
        pnodes[i].getNewCoord();

     maxerror = 0.0;
     for( size_t i = 0; i < nSize; i++)  {
         error = pnodes[i].setNewCoord();
         maxerror = max( error, maxerror );
     }
     if( maxerror < 1.0E-06) break;
     }
}

void CentroidalVoronoiTessalator :: smooth( Mesh *m, int iter )
{
     maxiter = iter;
     exactinit();

     mesh = m;
     int relexist0 = mesh->build_relations(0, 0);
     int relexist2 = mesh->build_relations(0, 2);

     mesh->search_boundary();

     SwapTriEdge edgeswap(mesh);

     assert( mesh->isHomogeneous() == 3 );

     mopt.shape_optimize(mesh);
     edgeswap.apply_rule( 1 );
     mopt.shape_optimize(mesh);

     for( int j = 0; j < 10; j++) {
         edgeswap.apply_rule();
     }

     if (!relexist0) mesh->clear_relations(0, 0);
     if (!relexist2) mesh->clear_relations(0, 2);
}

#endif

int main(int argc, char **argv )
{
     Jaal::Mesh *mesh = new Jaal::Mesh;

     assert( argc == 3 );
     mesh->readFromFile( argv[1] );

     MetisPartitioner mp;
     int npart = 10;
     mp.getPartition(mesh, npart);

     mesh->saveAs( argv[2] );
     exit(0);

     CentroidalVoronoiTessalator   cvt;
     cvt.smooth( mesh );

     mesh->saveAs( argv[2] );

     return 0;
}
