#pragma once


#include "Mesh.hpp"
#include "basic_math.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

class HeatConduction3D
{
    using Mat4D = Eigen::Matrix<double,4,4>;

public:
    HeatConduction3D();

    void setMesh( JMeshPtr m)
    {
        mesh = m;
        matrix_ready = 0;
    }

    void setBoundaryCondition( JNodePtr vtx, int, double val)
    {
        if( vtx == nullptr ) return;
        vtx->setAttribute("Dirichlet", val);
    }

    void setAttributeName( const std::string &n);
    void setMaxIterations( int n )
    {
        numIterations = n;
    }

    int solve();
private:
    JMeshPtr mesh;
    int    numBoundConditions;
    int    numIterations;
    bool   matrix_ready;
    std::string attribname;
    Eigen::SparseMatrix<double> M, A;
    Eigen::VectorXd  rhs, temperature;

    void build_global_matrix();
    void get_local_tet_mat(Vec4D &x, Vec4D &y, Vec4D &z, Mat4D &lmat);
    void apply_boundary_conditions();
    void solve_linear_system();
};

