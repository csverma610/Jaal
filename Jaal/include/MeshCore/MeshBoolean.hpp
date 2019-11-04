#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include <Eigen/Core>

#include <igl/boolean/mesh_boolean.h>

class JMeshBoolean
{
public:
    static const int MESH_UNION = 0;
    static const int MESH_COMPLEMENT   = 1;
    static const int MESH_INTERSECTION = 2;
    static const int MESH_DIFFERENCE   = 3;
    static const int MESH_SYMMETRIC_DIFFERENCE   = 4;
    static int  getOp( const string &op);

    void setMesh(const JMeshPtr &A, const JMeshPtr &b);
    JMeshPtr   getMesh(int op);
private:

    Eigen::MatrixXd VA, VB, VC;
    Eigen::MatrixXi FA, FB, FC;
};
