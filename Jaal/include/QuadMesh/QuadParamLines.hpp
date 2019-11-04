#pragma once

#include "Mesh.hpp"
#include "MeshTopology.hpp"

namespace Jaal
{
class  QuadParametricLines
{
public:
    QuadParametricLines()
    {
    }

    int setLines( JMeshPtr m );

private:
    void basicOp( JEdgePtr edge);

    std::deque<JEdgePtr> edgeQ;
    JMeshPtr mesh;
    JNodeSequence sequence;
};

}

