#pragma once

#include "Mesh.hpp"

using namespace Jaal;

class JMorseAnalysis
{
public:
    JMorseAnalysis()
    {
        mesh = NULL;
    }

    void setMesh(JMeshPtr m )
    {
        mesh = m;
    }

    void setHeightDirection( int dir);
    int  getCriticalNodes( JNodeSequence &n);

private:
    JMeshPtr mesh;
};

