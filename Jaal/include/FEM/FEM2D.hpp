#pragma once

#include <Eigen/Core>

#include "MeshCore/Mesh.hpp"
#include "MeshCore/SparseMatrix.hpp"
#include "MeshCore/StopWatch.hpp"
#include "FEM/FEMSpace.hpp"

/*
#include "MeshTangle.hpp"
#include "MeshExporter.hpp"
#include "BarycentricCoords.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
*/

int  getTexCoords1(const JFacePtr &face, const Point2D &xy, Point2D &uv);
int  getTexCoords2(const JFacePtr &face, const Point2D &xy, Point2D &uv);

bool isSymmetric( const Eigen::SparseMatrix<double> &A);
bool isSymmetric( const Eigen::MatrixXd &A);
void latexMatrix( const Eigen::MatrixXd &A);
void latexVector( const std::vector<double> &v);

////////////////////////////////////////////////////////////////////////////////

struct JFieldFunction {
    virtual double getScalar(const Point2D &) const = 0;
//   virtual Vec2D  getVector( const Point2D &p) const {}
};

////////////////////////////////////////////////////////////////////////////////

class JFEM2D
{
public :
    JFEM2D();

    void setMesh( const JMeshPtr &m);

    void setOrder( int e);

    void checkTangle( int v)
    {
        detectTangle = v;
    }

    void useAbsoluteJacobian( bool b )
    {
        absJacobian = b;
    }
    void setSolverMethod( int v);

    void setFieldFunction(const std::shared_ptr<JFieldFunction> &f)
    {
        fieldFunction = f;
    }

    virtual int solve();

    void    getXYCoords( const JFacePtr &face, const Point2D &uv, Point2D &xy);
    int     zeroJacobian( const JFacePtr &face, Point2D &pfirst, Point2D &psecond);
    double  getJacobian( const JFacePtr &f, const Point2D &uv);
    int     getSelfTangleRegion(const JFacePtr &f, vector<Point2D> &uvCoords, vector<Point2D> &xyCoords);

    int     splitSelfTangledQuad( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads);

protected:
    JMeshPtr mesh;

    int      shapeOrder;      //  linear or Quadrilateral)
    int      solverMethod;    // Direct, indirect etc...
    int      dofPerNode;      // 2 for 2D problems (u,v) in Ku = f matrix
    int      numDOF;

    bool     absJacobian;
    int      globalStatus;
    int      verbose;

    int      numTangledPairs;
    bool     detectTangle;
    int      numNodesPerElement[3];
    int      numGaussPoints[3];

    JStopWatch stopWatch;

    std::unique_ptr<JFEMSpace>  femSpace, tangleSpace;
    std::shared_ptr<JFieldFunction>  fieldFunction;

    Eigen::SparseMatrix<double> K, Ktangle;

    // Positions of Gauss points on canonical edge;
    std::vector<Point2D> edgeGaussPoints;
    Eigen::MatrixXd        Ne,  egradU;

    // Positions of Gauss points on canonical face;
    std::vector<Point3D> triGaussPoints;
    std::vector<Point3D> quadGaussPoints;

    std::vector<int>      freeDof;  // a map of fixed (u or v) values ...
    std::map<int,double>  dirichlet, neumann;  // which dofs are fixed

    virtual void applyDirichletBC() {}
    virtual void applyNeumannBC() {}
    virtual void applyBC() {}

    void getShapeFunc( const JFacePtr &face, const Point2D &uv,  std::vector<double> &N);
    void getShapeDeriv( const JFacePtr &face, const Point2D &uv, Eigen::MatrixXd &gradN);

    void getInvJacobian(const JFacePtr &face, Eigen::Matrix2d &invJ, double &dJ);
    void reduceSystem( Eigen::SparseMatrix<double> &A, std::vector<double> &b,
                       std::map<int,double> &fixMap, vector<int> &freeDof);

//    int   getJacobian( const JFacePtr &f, Matrix2d &mat, double &J);

    int   getNumOfNodes(const JEdgePtr &e) const;
    int   getNumOfNodes(const JFacePtr &e) const;

    int   initMesh();
    void  getCoords( const JFacePtr &face, Eigen::MatrixXd &xy);
    void  getNodes( const JFacePtr &face, std::vector<int> &nodes);
    void  mapUV(int triangle, const Point2D &uvtri, Point2D &uvquad);
    void  getAssemblyIndices(const JFacePtr &face, std::vector<int> &globalPos);
    void  getAssemblyIndices(const JEdgePtr &edge, std::vector<int> &globalPos);
    void  assembleK();
    void  assembleBC();
    void  genQuadraticNodes();
    int   getUVCoords( const JFacePtr &face, const Point2D &xy, Point2D &uv);
    bool  isTangled();
    void  detectTangleElements();
    void  zeroSearch( const JFacePtr &face, Point2D &p0, Point2D &p1, Point2D &pm, int &found);
    double getTriJacobian( const std::vector<Point2D> &triPnts);

//  Step 2 ...
    virtual int  buildLinearSystem() {}
    virtual int  solveLinearSystem() {}

    int splitSelfTangledCorner0( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads);
    int splitSelfTangledCorner1( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads);
    int splitSelfTangledCorner2( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads);
    int splitSelfTangledCorner3( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads);

};

