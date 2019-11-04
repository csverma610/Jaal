#pragma once

#include "Mesh.hpp"
#include "basic_math.hpp"

#ifdef USE_HASHMAP
#include <tr1/unordered_map>
#endif

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////

class  DijkstraShortestPath
{
public:

    void setMesh(JMeshPtr &m, const JMeshFilterPtr &f)
    {
        vsrc = 0;
        vdst = 0;
        mesh = m;
        filter = f;
        complete_path_known = 0;
    }

    JNodeSequence getPath( const JNodePtr &vs, const JNodePtr &vd);
    JNodeSequence getPath( const JNodePtr  &vs);

private:
    // Input parameters ...
    JMeshPtr mesh;
    JNodePtr vsrc, vdst;
    JMeshFilterPtr  filter;

    // Output data ...
    JNodeSequence  nodepath;

    // Local data ...
    int   distance_measure_method;
    bool  complete_path_known;

    struct LVertex {
        LVertex()
        {
            distance = 0.0;
            previous = nullptr;
            vertex   = nullptr;
        }

        size_t getID() const
        {
            return vertex->getID();
        }

        bool operator > ( const LVertex &rhs) const
        {
            return distance > rhs.distance;
        }

        bool operator < ( const LVertex &rhs) const
        {
            return distance < rhs.distance;
        }
        double   distance;   // Shortest distance from the source to this point
        JNodePtr vertex;    // Current Vertex
        JNodePtr previous;  // Previous vertex
    };

    /*
    #ifdef USE_HASHMAP
        std::tr1::unordered_map<JNodePtr, LVertex> vmap;
        std::tr1::unordered_map<JNodePtr,LVertex>::const_iterator miter;
    #endif
    */
    std::map<JNodePtr, LVertex> vmap;
    std::map<JNodePtr, LVertex>::const_iterator miter;

    std::priority_queue<LVertex, vector<LVertex>, greater<LVertex> > vertexQ;

    int   atomicOp( LVertex &node);
    void  fastmarching();  // Fast Marching Style algorithm O(nlogn) with heap
    void  traceback();
    double getCost(const LVertex &vi, const LVertex &vj ) const;
};

