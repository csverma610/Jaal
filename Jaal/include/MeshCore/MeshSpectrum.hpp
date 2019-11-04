#pragma once

#include "Mesh.hpp"
#include <Eigen/Core>

/*
#include "MeshMatrix.hpp"
#include <armadillo>
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/repdiag.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
#include <iostream>
*/
class JMeshSpectrum
{
    typedef std::tuple<int,int, double> MatrixData;
public:
    void setMesh(const JMeshPtr &m)
    {
        mesh = m;
    }

    void genEigenVectors(int n);

    std::vector<double> getEigenVector(int i) const;
private:
    JMeshPtr mesh;
    Eigen::MatrixXd EVec;
};
