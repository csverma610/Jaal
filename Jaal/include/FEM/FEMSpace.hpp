#pragma once

#include "MeshCore/Mesh.hpp"

struct JFiniteElement {
    JFiniteElement() {}

    void setOrder(int n)
    {
        order = n;
    }
    int  order = 1;
    virtual int getNumNodes() const = 0;
    virtual int getNumGaussPoints() const = 0;
};

////////////////////////////////////////////////////////////////////////////////

struct JEdgeElement : public JFiniteElement {
    JEdgeElement(int n = 2)
    {
        setGaussPoints(n);
    }

    int   getNumNodes() const
    {
        if( order == 1) return 2;
        if( order == 2) return 3;
        return 0;
    }
    int   getNumGaussPoints() const
    {
        return gaussPoints.size();
    }

    int  setGaussPoints(int n);

    const vector<Point2D> &getGaussUPoints() const
    {
        return gaussPoints;
    }

    int  getShapeFunc( double u,   vector<double> &N);
    int  getShapeDeriv( double u,  vector<double> &gradN);

    double  eval(double xi, const vector<double> &phi);
    Array2D eval(double xi, const vector<Array2D> &phi);
    Array3D eval(double xi, const vector<Array3D> &phi);

    void  getGaussXY(const vector<Point2D> &xyPoints, vector<Point2D> &xygauss);

private:
    vector<Point2D> gaussPoints;
};
////////////////////////////////////////////////////////////////////////////////

struct JTriElement : public JFiniteElement {
    JTriElement(int n = 3)
    {
        setGaussPoints(n);
    }

    int   getNumNodes() const
    {
        if( order == 1) return 3;
        if( order == 2) return 6;
    }

    int   getNumGaussPoints() const
    {
        return gaussPoints.size();
    }

    int    setGaussPoints(int n);

    const vector<Point3D> &getGaussUVPoints() const
    {
        return gaussPoints;
    }

    int   getShapeFunc( const Point2D &uv, vector<double> &N);
    int   getShapeDeriv( const Point2D &uv,  Eigen::MatrixXd &gradN);
    double getJacobian( const vector<Point2D> &xy, const Point2D &uv);

    double  eval(const Point2D &uv, const vector<double> &phi);
    Array2D eval(const Point2D &uv, const vector<Array2D> &phi);
    Array3D eval(const Point2D &uv, const vector<Array3D> &phi);

    int  getGaussXYPoints(const vector<Point2D> &xyPoints, vector<Point2D> &xygauss);
    int  getXYCoords(const vector<Point2D> &xyPoints,  const Point2D &uv, Point2D &xy);
    int  getUVCoords(const vector<Point2D> &xyPoints, const Point2D &xy, Point2D &uv);

private:
    vector<Point3D> gaussPoints;
    int    getT3ShapeFunc(  const Point2D &uv, vector<double> &N);
    int    getT6ShapeFunc(  const Point2D &uv, vector<double> &N);
    int    getT3ShapeDeriv( const Point2D &uv, Eigen::MatrixXd  &gradN);
    double getT3Jacobian( const vector<Point2D> &p, const Point2D &uv);
};

////////////////////////////////////////////////////////////////////////////////

struct JQuadElement : public JFiniteElement {
    JQuadElement(int n = 4)
    {
        setGaussPoints(n);
    }

    int   getNumNodes() const
    {
        if( order == 1) return 4;
        if( order == 2) return 8;
    }

    int   getNumGaussPoints() const
    {
        return gaussPoints.size();
    }

    int  setGaussPoints(int n);

    const vector<Point3D> &getGaussUVPoints() const
    {
        return gaussPoints;
    }

    int  getShapeFunc( const Point2D &uv, vector<double> &N);
    int  getGaussXYCoords( const vector<Point2D> &xyCoord, vector<Point2D> &p, int n = 0);

    int  getShapeDeriv( const Point2D &uv,  Eigen::MatrixXd &gradN);

    int  getXYCoords( const vector<Point2D> &xyCoords, const Point2D &uv, Point2D &xy);
    int  getUVCoords( const vector<Point2D> &xyCoords, const Point2D &xy, Point2D &uv);
    double getJacobian( const vector<Point2D> &xy, const Point2D &uv);

    double  eval(const Point2D &uv, const vector<double> &phi);
    Array2D eval(const Point2D &uv, const vector<Array2D> &phi);
    Array3D eval(const Point2D &uv, const vector<Array3D> &phi);

    void  getXYPoint(const vector<Point2D> &xyPoints,  const Point2D &uv, Point2D &xy);
    void  getUVPoints(const vector<Point2D> &xyPoints, const Point2D &xy, Point2D &uv);

private:
    vector<Point3D> gaussPoints;
};

////////////////////////////////////////////////////////////////////////////////


struct JConcaveQuadElement {
    int  getGaussUVCoords( const vector<Point2D> &xyCoords, vector<Point3D> &p, int n = 4);
    int  getGaussXYCoords( const vector<Point2D> &xyCoords, vector<Point3D> &p, int n = 4);

    int  getShapeFunc1( int corner, const Point2D &ab, const Point2D &uv, vector<double> &N);
    int  getShapeDeriv1(int corner, const Point2D &ab, const Point2D &uv, Eigen::MatrixXd  &gradN);

    double getJacobian( const vector<Point2D> &p, const Point2D &uv);

    // For concave quad element, there may be positive and negative jacobians. Identify
    // the line of zero Jacobian..
    int  getZeroJacobianLine( const vector<Point2D> &poly, Point2D &p1, Point2D &p2);
    int  getReflectionPoint( const vector<Point2D> &xyCoords, Point2D &ab);

    int getPositiveJacobianRegion( const vector<Point2D> &poly, vector<Point2D> &region);
    int getNegativeJacobianRegion( const vector<Point2D> &poly, vector<Point2D> &region);
    int getConcaveCorner( const vector<Point2D> &poly);

    int split( const vector<Point2D> &poly, vector<Point2D> &uv, vector<Array4I> &quads);

private:
    JQuadElement stdQuad;

    void zeroSearch( const vector<Point2D> &points, Point2D &p0, Point2D &p1, Point2D &pm, int &found);

    int getReflection0( const vector<Point2D> &poly, Point2D &uv);
    int getReflection1( const vector<Point2D> &poly, Point2D &uv);
    int getReflection2( const vector<Point2D> &poly, Point2D &uv);
    int getReflection3( const vector<Point2D> &poly, Point2D &uv);

    int splitCorner0( const vector<Point2D> &poly, vector<Point2D> &uv, vector<Array4I> &quads);
    int splitCorner1( const vector<Point2D> &poly, vector<Point2D> &uv, vector<Array4I> &quads);
    int splitCorner2( const vector<Point2D> &poly, vector<Point2D> &uv, vector<Array4I> &quads);
    int splitCorner3( const vector<Point2D> &poly, vector<Point2D> &uv, vector<Array4I> &quads);

    int getCorner0Func( const Point2D &ab, const Point2D &uv, vector<double> &N);
    int getCorner1Func( const Point2D &ab, const Point2D &uv, vector<double> &N);
    int getCorner2Func( const Point2D &ab, const Point2D &uv, vector<double> &N);
    int getCorner3Func( const Point2D &ab, const Point2D &uv, vector<double> &N);

    int getCorner0Deriv( const Point2D &ab, const Point2D &uv, Eigen::MatrixXd &gradN);
    int getCorner1Deriv( const Point2D &ab, const Point2D &uv, Eigen::MatrixXd &gradN);
    int getCorner2Deriv( const Point2D &ab, const Point2D &uv, Eigen::MatrixXd &gradN);
    int getCorner3Deriv( const Point2D &ab, const Point2D &uv, Eigen::MatrixXd &gradN);
};

////////////////////////////////////////////////////////////////////////////////

class JFEMSpace
{
public:
    static const int  EDGE         = 1;
    static const int  TRIANGLE     = 23;
    static const int  QUAD         = 24;

    static const int  TETRAHEDRON  = 34;
    static const int  HEXAHEDRON   = 38;

    JFEMSpace();

    void setOrder( int n);

    JEdgeElement  edgeElement;
    JTriElement   triElement;
    JQuadElement  quadElement;
};

///////////////////////////////////////////////////////////////////////////////
int JQuadBaryShapeFunc( const vector<Point2D> &quadPoints, const Point2D &uv, vector<double> &N);
