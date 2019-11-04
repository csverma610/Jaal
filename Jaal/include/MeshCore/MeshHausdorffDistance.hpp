#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include "AllTetMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"
#include <igl/point_mesh_squared_distance.cpp>

class JMeshHausdorffDistance
{
public:
    JMeshHausdorffDistance() {}

    void setSource( const JMeshPtr &m)
    {
        srcMesh = m;
    }
    void setTarget( const JMeshPtr &m)
    {
        dstMesh = m;
    }

    vector<double>  getDistance(int dir);

private:
    JMeshPtr  srcMesh, dstMesh;
};
