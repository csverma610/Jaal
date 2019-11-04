#pragma once

#include "Mesh.hpp"

#include <Eigen/Core>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/lscm.h>
#include <igl/harmonic.h>
#include <igl/arap.h>
#include "MeshMatrix.hpp"
#include "MeshAffineTransforms.hpp"

//#include <igl/svd3x3/arap.h>

// Code borrowed from : http://www.riken.jp/brict/Yoshizawa/Research/Param.html

class JSurfaceParameterization
{
public:
    // Weighttype ...
    static const int SHAPE_PRESERVING = 0;
    static const int TUTTE            = 1;
    static const int HARMONIC_MAP     = 2;
    static const int INTRINSIC_MAP    = 3;
    static const int MEAN_VALUE       = 4;
    static const int LEAST_SQUARE_CONFORMAL   = 5;
    static const int AS_RIGID_AS_POSSIBLE  = 6;

    // BoundaryType ...
    static const int SQUARE           = 0;
    static const int CIRCLE           = 1;
    static const int NATURAL_BOUNDARY = 2;

    JSurfaceParameterization();

    int  setPatchMesh( const JMeshPtr &m);

    void setWeight( int w)
    {
        weighttype = w;
    }
    void setBoundary(int b)
    {
        boundarytype = b;
    }
    void setStartCornerID( int id)
    {
        vertexCornerID = id;
    }
    void setGamma(double g)
    {
        gammaP = g;
    }
    void setNumIterations( int n)
    {
        numIters = n;
    }

    JMeshPtr  getParamMesh();

private:
    JMeshPtr  patchMesh;
    int    weighttype;
    int    boundarytype;
    int    numIters;
    int    smooth;
    int    vertexCornerID;
    double gammaP;
    double intrinsiclambda;

    bool isDisk();
    int  writePly2File();
    void normalize( const JMeshPtr &m);
    JMeshPtr readPly2File();
    JMeshPtr LeastSquareConformalMap();
    JMeshPtr AsRigidAsPossible();
};

