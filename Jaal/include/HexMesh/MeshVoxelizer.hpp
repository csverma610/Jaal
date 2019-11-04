#pragma once

#include "Mesh.hpp"
#include "AllHexMeshGenerator.hpp"

class JMeshVoxelizer
{
public:
    JVoxelMeshPtr  genVoxels(const JMeshPtr &mesh, int size = 128);
    JVoxelMeshPtr  genVoxels(const JNodeSequence &nodes, int nsize = 128);

private:
    JMeshPtr       bgmesh;
    JVoxelMeshPtr  voxmesh;

    double length[3],  origin[3], gridSpacing[3];
    int    nodeDim[3], cellDim[3];

    void checkOverlap(const JFacePtr &face);
    void checkOverlap(double triVerts[3][3], int i, int j, int k);
};

