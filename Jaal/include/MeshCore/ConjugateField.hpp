#pragma once

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/conjugate_frame_fields.h>
#include <igl/ConjugateFFSolverData.h>
#include <igl/dot_row.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/n_polyvector.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <vector>
#include <cstdlib>

#include "SurfaceVectorField.hpp"

class JConjugateField : public JSurfaceVectorField
{
public:
    void setMesh( const JMeshPtr &m);
    int  genField();
    void genRandomConstraints(int nRandom);

private:
    Eigen::MatrixXd smooth_pvf;
    Eigen::MatrixXd conjugate_pvf;
    Eigen::VectorXd conjugacy_s;
    Eigen::VectorXd conjugacy_c;
    igl::ConjugateFFSolverData<Eigen::MatrixXd, Eigen::MatrixXi> *csdata;
};
