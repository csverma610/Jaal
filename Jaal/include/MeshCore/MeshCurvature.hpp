#pragma once

#include "Mesh.hpp"
#include <Eigen/Core>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/principal_curvature.h>
#include "MeshMatrix.hpp"

class JMeshCurvature
{
public:
    void setMesh( const JMeshPtr &m);

    void setEvalPos( int p)
    {
        evalPos = p;
    }

    int  setGaussianCurvature();
    int  setMeanCurvature(int method = 0);
    void setScale( double s)
    {
        scale = s;
    }

    std::pair<JMeshPtr, JMeshPtr> getCurvatureDirections();

private:
    JMeshPtr mesh;
    double   scale;
    int      evalPos;
    JMeshPtr getVectors( const Eigen::MatrixXd &headPoints,
                         const Eigen::MatrixXd &tailPoints);
    void   assignAverage( const JFacePtr &f,  const string &s);
};

