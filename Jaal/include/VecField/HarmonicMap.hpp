#pragma once

#include "MeshInterpolation.hpp"

#include "Mesh.hpp"
#include "AllTetMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"
#include <igl/harmonic.h>

class JHarmonicMap
{
public:
    void setSource( const JMeshPtr &m);
    void setTarget( const JMeshPtr &m);
    JMeshPtr getDeformedMesh( double t = 1.0);

private:
    int order = 2;
    JMeshPtr simplicialMesh;
    JMeshPtr deformedMesh;

    Eigen::MatrixXd V, bc;
    Eigen::MatrixXi F;
    Eigen::VectorXi b;

    void init2D( const JMeshPtr &m);
    void init3D( const JMeshPtr &m);
    void solveSystem();
    void setBoundConditions( double t);
};
