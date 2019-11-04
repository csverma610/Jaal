#pragma once

#define EIGNN_SUPERLU_SUPPORT

#include <map>

#include "FEMSpace.hpp"
#include "MeshTangle.hpp"
#include "LocallyInjectiveMap.hpp"
#include "MeshExporter.hpp"
#include "SparseMatrix.hpp"
#include "StopWatch.hpp"
#include "BarycentricCoords.hpp"
#include "TangleElements2D.hpp"
#include "FEM2D.hpp"

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

int getConcaveCorner( const JFacePtr &f);

class JElasticity2D : public JFEM2D
{
public:
    static const int  PLAIN_STRESS = 0;
    static const int  PLAIN_STRAIN = 1;

    JElasticity2D() ;

    void setYoungModulus( double  v)
    {
        E  = v;
        computeDMatrix();
    }

    void setPoissonRatio( double  v)
    {
        nu = v;
        computeDMatrix();
    }

    void setXForce( const JEdgeSequence &edges, double force);
    void setYForce( const JEdgeSequence &edges, double force);

    void fixX( const JNodePtr &vtx, double u);
    void fixY( const JNodePtr &vtx, double v);
    void fixX( const JEdgeSequence &boundaryEdges, double u = 0.0);
    void fixY( const JEdgeSequence &boundaryEdges, double v = 0.0);
    void fixEdges( const JEdgeSequence &boundaryEdges, double uv = 0.0);

    void setForceVector( double *v)
    {
        F[0] = v[0];
        F[1] = v[1];
    }

    void setBodyForce( double  *f)
    {
        bodyForce[0] = f[0];
        bodyForce[1] = f[1];
    }

    void saveAs( const string &filename, const string &varname);
    void getDisplacements( vector<double> &u, vector<double> &v);

    // From the displacements, calculate the stresses ....
    void computeStress();

    int getValueAt( char c, int id, double &val)
    {
        if( c == 'U') {
            val = u[id];
            return 0;
        }

        if( c == 'V') {
            val = v[id];
            return 0;
        }

        if( c == 'S') {
            val = nodeVonMises[id];
            return 0;
        }
        return 1;
    }

private:
    int      problem;
    double   E;               // Young modulus;
    double   nu;              // Poission Ratio.
    double   F[2], bodyForce[2];

    Matrix3d  D;
    vector<double> u, v;  // Displacement vector in Ku = f
    vector<double> f, freduced;  // force vector in Ku = f;

    vector<vector<double> > NFace; // Basis function for a face.
    vector<MatrixXd> gradUV;     // Gradident of Basis function for each Gauss point..
    vector<double>   faceArea;     // Area of each face ( > 0 )
    vector<double>   nodeVonMises, faceVonMises;
    vector<Matrix2d> faceStrain, faceStress;

    void initParams();

    // Given E and nu, calculate D matrix once ...
    void  computeDMatrix();

    void BMatrix( const JFacePtr &face, const Point2D &uv, MatrixXd &B);
    void BMatrix( const JFacePtr &face, const Point2D &uv, MatrixXd &B, double &dJ);

    void BMatrix(const double *uv, const Matrix2d &invJ, MatrixXd &B);
    void BMatrix(const MatrixXd &gradUV, const Matrix2d &invJ, MatrixXd &Bg);

    // Assemble the global matrix from each element for Ku = f.
    void  assembleK();

    void applyNeumannBC();
    void applyDirichletBC();
    void applyBC();

    // Calculate the right hand side and modify the Left hand side K Matrix...
    void  assembleBC();

    //   Integration over the element for K = sum( B^T*D*B dx)
    int   integrate(const JFacePtr &face, const JFacePtr &f, MatrixXd &K);
    void  integrate(const JFacePtr &face, MatrixXd &K, vector<double> &f);
    void  convexIntegrate(const JFacePtr &face, MatrixXd &K, vector<double> &f);
    void  concaveIntegrate(const JFacePtr &face, MatrixXd &K, vector<double> &f);
    void  meanValueIntegrate(const JFacePtr &face, MatrixXd &K, vector<double> &f);
    void  concaveQuadTriIntegrate(const JFacePtr &quadface, int tid, MatrixXd &Ke, vector<double> &fe);
    void  meanValueTriIntegrate(const JFacePtr &quadface, int tid, MatrixXd &Ke, vector<double> &fe);
    void  integrate(const JEdgePtr &e, double val, vector<double> &f);
    void  integrate( const JFacePtr &face1, const JFacePtr &face2,
                     const vector<Point2D> &triPnts, MatrixXd &Ke);

    void computeStress( const JFacePtr &f);

    void assembleTangledK( const JFacePtr &f0, const JFacePtr &f1);
    void assembleTangledK();

    int  buildLinearSystem();
    int  solveLinearSystem();

};

//////////////////////////////////////////////////////////////////////////////

