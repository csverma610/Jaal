#pragma once

#include "Mesh.hpp"
#include <ANN/ANN.h>

class JNearestNeighbours
{
public:
    JNearestNeighbours();
    ~JNearestNeighbours();

    void setMesh( const JMeshPtr &m);
    void setCloud( const JNodeSequence &seq)
    {
        inCloud = seq;
    }

    JNodeSequence  nearestKSearch( const Point3D &p, int k);
    JNodePtr      nearest( const Point3D &p);

    JNodeSequence  radiusSearch( const Point3D &p, int k);

private:
    JMeshPtr mesh;
    JNodeSequence  inCloud;

    ANNpointArray     dataPts;
    ANNpoint          queryPoint;
    ANNidxArray       nnIdx;
    ANNdistArray      dists;
    std::scoped_ptr<ANNkd_tree> kdTree;

    void buildTree();
};
