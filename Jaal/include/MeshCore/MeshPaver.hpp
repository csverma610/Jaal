#pragma once

#include "Mesh.hpp"
#include "MeshTopology.hpp"

using namespace Jaal;

class JMeshPaver
{
public:
    JMeshPaver()
    {
        hybrid = 0;
    }

    void setMesh(JMeshPtr m)
    {
        mesh = m;
    }
    void setNumBoundaryLayers( int n)
    {
        numLayers = n;
    }
    void setHybrid( bool v)
    {
        hybrid = v;
    }

    void execute();

private:
    JMeshPtr mesh;
    int    numLayers;
    bool   hybrid;
    int    elemType;

    JNodeSequence boundnodes;
    JEdgeSequence boundedges;

    void   init();
    void   addNewLayer();
    void   optimize();
};
