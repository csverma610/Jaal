#pragma once

#include "MeshCore/Mesh.hpp"

class  JQuadExtractor
{
public:
    JQuadExtractor()
    {
        uvscale = 1.0;
    }
    void setMesh( const JMeshPtr &m)
    {
        mesh = m;
    }
    void setUVScale( int n)
    {
        uvscale = n;
    }
    size_t getNumIntegerUV();
    JMeshPtr  getQuadMesh();
private:
    JMeshPtr mesh;
    int uvscale;
};

class AllQuadMeshGenerator
{
public:
    static const int   GREEDY_MATCHING = 0;
    static const int   EDMONDS_MATCHING = 1;
    static const int   BINARY_TREE_MATCHING = 2;

    static const int   TRI2QUAD_BASE_QUADS  = 0;
    static const int   TRIMATCH_BASE_QUADS  = 1;
    static const int   SKELETON_BASE_QUADS  = 2;
    static const int   CONVEX_BASE_POLYGONS = 3;
    static const int   BIPARTITE_BASE_POLYGONS = 4;

    static JMeshPtr SchneiderPyramid();
    static JMeshPtr getSquareHoles(int nx= 1, int ny = 1, double innerRadius = 1, double oRadius = 2);
    static JMeshPtr getStructuredMesh(int *dim, double *len = nullptr, double *org= nullptr, bool texCoord = 0);
    static JMeshPtr getSierpinski(int nlevel );

    void   setMesh(const JMeshPtr &m);

    JMeshPtr refineAll();
    JMeshPtr minimalRefine();
    JMeshPtr getSimpleTris2Quads();
    JMeshPtr getHamiltonianQuads();
    JMeshPtr getTrianglesMatching(int algo = GREEDY_MATCHING);
    JMeshPtr getSimpleQuadMesh();

    JMeshPtr getBaseQuadMesh(int algo);

    int  makeEvenTriangulation();
    int  splitBoundQuad2QuadTriangle(const JEdgePtr &e);

private:
    JMeshPtr  mesh, trimesh, quadmesh;
    int   getEdmondsGraphMatching();
    int   getBinaryTreeMatching();
    int   getGreedyTrianglesMatching();

    int refineTriangle( const JFacePtr &f);
    int refineQuad( const JFacePtr &f);
    int verify_matching(const JMeshPtr &m, const vector<JFacePair> &matching);
    JMeshPtr collapse_matched_triangles(const JMeshPtr &mesh, const vector<JFacePair> &matching);

    vector<Point3D> midPoints;
    JNodeSequence skippedNodes;
    JMeshPtr getBaseTriMesh( bool resample );

    JMeshPtr getBaseQuads_Tri2Quads();
    JMeshPtr getBaseQuads_TriMatch();
    JMeshPtr getBaseQuads_Skeleton();
    JMeshPtr getBaseQuads_EdgeMatch();

    JMeshPtr getBasePolygons();

};
