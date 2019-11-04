#pragma once

#include "Mesh.hpp"

using namespace Jaal;

class JBoySurface
{
public:
    JMeshPtr getApery(int nx, int ny);
    JMeshPtr getBryant(int nx, int ny);
    JMeshPtr getMin1();
    JMeshPtr getMin2();
    JMeshPtr getMin3();
private:
    JMeshPtr newmesh;
    JNodeSequence nodes;
    void createNodes();
    void getQuadMesh(int nx, int ny);
};


