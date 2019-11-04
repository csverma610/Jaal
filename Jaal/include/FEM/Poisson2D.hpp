#pragma once

#define EIGNN_SUPERLU_SUPPORT

#include "FEM2D.hpp"

#include <Eigen/Core>

/*
#include <map>
#include "FEMSpace.hpp"
#include "MeshTangle.hpp"

#include "MeshExporter.hpp"
#include "SparseMatrix.hpp"
*/

///////////////////////////////////////////////////////////////////////////////
//
// Description:  For the 2D triangle mesh, this class calculates the displacement,
//               and stresses.
// What input is needed:
//	(A)  A triangle mesh and mesh need not be simplicial in the sense that
//           triangle elements can overlap. But there are some conditiosns:
//	     (1) The mesh shhould be sliver free.
//           (2) No element should tangle with boundary edges.
//      (B)  Boundary conditions:
//           Force boundary conditions: Force applied to the geometric edges.
//           Displacemennt specfied to some to the boundary nodes. Most of the
//           the values are zero.
//
// License:  It is free and has absolutely no restrictions whatsoever.
///////////////////////////////////////////////////////////////////////////////
//
// The initial code was written in Matlab, but was converted to C++ because of
// my inertia to learn new language and (II) I want to experiment with new C++
// features.
//          Chaman Singh Verma
//          Department of computer Sciences.
//          University of Wisconsin,
//          Madison
// The original Matlab code was provided by:
//          Prof. Suresh Krishnan,
//          Department of Mechanical Engineering.
//          University of Wisconsin, Madison
//          Madison.
//
//////////////////////////////////////////////////////////////////////////////

class JPoisson2D : public JFEM2D
{
public:
    JPoisson2D() ;

    void setDirichletValue( const JNodePtr &v, double u);
    void setDirichletValue( const JEdgePtr &e, double u);
    void setDirichletValue( const JEdgeSequence &ve, double u);

    void setNeumannValue( const JNodePtr &v, double u);
    void setNeumannValue( const JEdgePtr &e, double u);
    void setNeumannValue( const JEdgeSequence &ve, double u);

    void getField( vector<double> &u);

    int getValueAt( int id, double &val)
    {
        val = u[id];
        return 0;
    }

    int solve();


    void saveAs( const string &filename, const string &varname);

    void  integrate_self_tangle(const JFacePtr &face, Eigen::MatrixXd &Ke, std::vector<double> &fe);

private:
//  JMeshTangle meshtangle;
//  vector< vector<double> >  NEdge, gradU; //  Basis function and gradient at Gauss points for an edge.

    std::vector<double> u;             // Displacement vector in Ku = f
    std::vector<double> f, freduced;  // force vector in Ku = f;
    std::vector<std::vector<double> > NFace; // Basis function for a face.
    std::vector<Eigen::MatrixXd> gradUV;     // Gradident of Basis function for each Gauss point..

    void initParams();

    void BMatrix( const JFacePtr &face, const Point2D &uv, Eigen::MatrixXd &B);
    void BMatrix( const JFacePtr &face, const Point2D &uv, Eigen::MatrixXd &B, double &dJ);
    void BMatrix(const double *uv, const Eigen::Matrix2d &invJ, Eigen::MatrixXd &B);
    void BMatrix(const Eigen::MatrixXd &gradUV, const Eigen::Matrix2d &invJ, Eigen::MatrixXd &Bg);

    int   buildLinearSystem();
    void  assembleK();
    void  applyNeumannBC();
    void  applyDirichletBC();
    void  applyBC();

    // Integration of a triangles within the overlap region ...
    int   integrate(const JFacePtr &face, const JFacePtr &f, Eigen::MatrixXd &K);

    void  integrate( const JFacePtr &face1, const JFacePtr &face2,
                     const std::vector<Point2D> &triPnts, Eigen::MatrixXd &Ke);

    //   Integration over the element for K = sum( B^T*D*B dx)
    void  integrate(const JFacePtr &face, Eigen::MatrixXd &K, std::vector<double> &f);

    //  Integrate over the boundry edge for the right hand side
    //  vector f = sum ( N f dx )
    void  integrate(const JEdgePtr &e, double val, std::vector<double> &f);

    int  solveLinearSystem();

    void assembleTangledK( const JFacePtr &f0, const JFacePtr &f1);
    void assembleTangledK();
};

//////////////////////////////////////////////////////////////////////////////

