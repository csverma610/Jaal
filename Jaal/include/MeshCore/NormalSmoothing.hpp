#pragma once

#include "Mesh.hpp"

class JMeshNormalSmoothing
{
public:
    JMeshNormalSmoothing()
    {
        numIterations = 5;
    }

    void setMesh( const JMeshPtr &m)
    {
        mesh = m;
    }
    void setNumIterations(int n )
    {
        numIterations = n;
    }
    void execute();

private:
    JMeshPtr mesh;
    int numIterations;
    void atomicOp( const JNodePtr &p);
    void atomicOp( const JFacePtr &p);
};

