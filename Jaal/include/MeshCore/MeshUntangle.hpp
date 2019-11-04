#pragma once

#include "Mesh.hpp"
#include "MeshOptimization.hpp"
#include "MeshLaplacian.hpp"
#include "LocallyInjectiveMap.hpp"
#include "LloydOptimizer.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AllTetMeshGenerator.hpp"

class JMeshUntangle
{
public:
    JMeshUntangle();

    ~JMeshUntangle();

    // Input mesh containing tangled elements...
    void setMesh( const JMeshPtr &m);

    // How many elements are inverted in the tri/tet mesh ..
    size_t countInverted( const JMeshPtr &m);
    // Back projection energy ....
    void   setEnergyType( int e)
    {
        energyType = e;
    }

    // settting minimum number of inflation iterations ?
    void   setNumBoundaryModifications( int n)
    {
        maxBoundaryModifications = n;
    }
    void   setNumBackProjections( int n)
    {
        maxBackProjections = n;
    }

    // How many nodes on the boundary will be convexified in order
    // to make mesh tangle-free ?
    size_t getNumNodesConvexified() const;
    bool   setInflateAll( bool b)
    {
        inflateAll = b;
    }

    // Return the inflated mesh, which sould be tangle free ...
    JMeshPtr getInflatedMesh();
    JMeshPtr getSimplicialMesh();

    // Inflated mesh to be bring back to the original shape along with the
    // uninverted elements...
    void  backProject();

    // Perform inflate and projection in one call.
    void smoothAll()
    {
        getInflatedMesh();
        backProject();
    }

    JFaceSequence getInvertedFaces();
    JCellSequence getInvertedCells();

    vector<double>  getOffset() const
    {
        return offset;
    }

    double  getAreaChanged() const
    {
        return 100*(area[1] - area[0] )/area[0];
    }
    double  getVolumeChanged() const;
    double  getMaxDistance() const;

private:
    JMeshPtr  srcMesh, inflatedMesh, triMesh, tetMesh;
    vector<double> offset;
    double area[2], volume[2];

    int  entityDim;
    int  energyType;
    int  maxBoundaryModifications   = 100;
    int  maxBackProjections = 1;

    bool inflateAll    = 0;

    JLaplaceMeshSmoother  lap;
    JLloydMeshOptimizer   lloyd;
    JLocallyInjectiveMap  lim;
    JMeshNonlinearOptimization mopt;

    JNodeSet concaveCorners;
    int smoothCorner( const JNodePtr &v);
    int smoothBoundary();

    JMeshPtr getInflated2D();
    JMeshPtr getInflated3D();

    int  updateSource();
    int  clear();
};
