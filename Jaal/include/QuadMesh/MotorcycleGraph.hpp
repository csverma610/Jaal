#pragma once

#ifndef MOTORGRAPH_H
#define MOTORGRAPH_H

#include "Mesh.hpp"
#include "MeshTopology.hpp"
#include "MeshPartitioner.hpp"

namespace Jaal
{

struct JMotorPath {
    JNodeSequence nodes;
    int   getNumJunctions();
};


class JSingularityGraph
{
public:
    void setMesh( const JMeshPtr &m);

    JEdgeSequence getPaths( const JNodeSequence &v );
    JEdgeSequence getAllPaths();
    JMeshPtr     getGraph();

protected:
    static JLogger *logger;
    JMeshPtr mesh;
    JNodeSequence singularNodes;
    JNodeSequence nodeSeq;
    JEdgeSequence edgeSeq;
    JFaceSequence faceSeq;
    bool diskTopology;
    int  advance_single_quad_step();
    int  getBranches(const JNodePtr &v);
    int  create_new_branch( const JNodePtr &v0, const JNodePtr &v1);
    void getQuadPaths();
    void getQuadPaths( const JNodeSequence &v);
    void getTriPaths();
    void getTriPaths( const JNodeSequence &v);
    void advance_tri_path(const JNodePtr &v1, const JNodePtr &v2);
    int  makeCyclicNodes( const JNodePtr &v);
};


class  JMotorcycleGraph : public JSingularityGraph
{
public:
    const static int END_AT_TERMINALS = 0;
    const static int END_AT_CROSSINGS = 1;

    JMotorcycleGraph()
    {
        mesh = nullptr;
    }

    /////////////////////////////////////////////////////////////////////////////
    // There are two ways to advance along the track.
    // (1) Incremental:  All track propagate simulateneously and
    //                   at the intersection, only one is allowed
    //                   to proceed towards other irregular node.
    // (11) Greedy   :   Start from one vertex and complete its
    //                   track.
    // It is not clear which method is better, but "greedy" may likely
    // give very high aspect ratio quad patches. Incremental on the
    // other hand may produce many small patches..
    //
    /////////////////////////////////////////////////////////////////////////////

    JNodeSequence getJunctions();

    int  getPartitions();
    int  getPartitions( const JNodeSequence &v);

    JMeshPtr  getToplogy();

    int optimize();
private:

    int  graph_topology;

};

}

#endif

