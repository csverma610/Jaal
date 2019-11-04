#pragma once

#include "Mesh.hpp"

struct AllTriMeshGenerator {

    static JMeshPtr  getStructuredMesh(int nx, int ny);
    static JMeshPtr  getCylinder( const Point3D &p0, const Point3D &p1, double radius, int nR);
    static JMeshPtr  fromQuad2Tri(int *dim, double *len = nullptr, double *org= nullptr, bool texCoord = 0);
    static JMeshPtr  fromQuad4Tri(int *dim, double *len = nullptr, double *org= nullptr, bool texCoord = 0);
    static JMeshPtr  fromQuadMesh(const JMeshPtr &quadmesh, JNodeSequence &steiner, int type = 2); // or type == 4
    static JMeshPtr  fromQuadMesh(const JMeshPtr &quadmesh, int type = 2);
    static JMeshPtr  fromPolyMesh(const JMeshPtr &polymesh, JNodeSequence &steiner);
    static JMeshPtr  fromPolyMesh(const JMeshPtr &polymesh);
    static JMeshPtr  getSierpinski(int nlevel ) ;
    static int       getIsotropicMesh( const JMeshPtr &m) ;

    static void  getCongruentMesh(const JMeshPtr &mesh);
};

