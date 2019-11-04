#pragma once

#include "Mesh.hpp"
#include "NearestNeighbours.hpp"
#include "PolygonQuadMesher.hpp"
#include "MeshRefine.hpp"
#include "Doublet.hpp"

class JCirclePatchQuadMesh
{
public:
    void setMesh(const JMeshPtr &m);

    bool isEmpty()
    {
        return patchFaces.empty();
    }

    void setCircle( const JCircle &m);
    void setCenter( const Point3D &c);
    JCircle  getCircle() const
    {
        return circle;
    }

    int  getNumSingularities() const
    {
        return numSingularities;
    }
    int  getNumBoundaryNodes() const
    {
        return patchEdges.size();
    }

    JFaceSequence getFaces()
    {
        return patchFaces;
    }
    JEdgeSequence getBoundaryEdges()
    {
        return patchEdges;
    }
    JNodeSequence getBoundaryNodes()
    {
        return patchNodes;
    }

    void setAdaptationFactor( double a)
    {
        adaptFactor = a;
    }
    void setExpectedEdgeLength( double a)
    {
        expectedEdgeLength = a;
    }

    void build();
    void remesh();
    void clear();

    JFaceSequence getNewFaces() const
    {
        return newFaces;
    }

private:
    JMeshPtr mesh;

    JCircle  circle;

    JNodeSequence patchNodes, newNodes;
    JEdgeSequence patchEdges, newEdges;
    JFaceSequence patchFaces, newFaces;

    double meanArea;
    int    numSingularities = 0;
    double adaptFactor = 1.0;
    double expectedEdgeLength = 0.0;

    std::unique_ptr<JNearestNeighbours> nearSearch;
};
