#include "AllQuadMeshGenerator.hpp"
#include "StopWatch.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include "EdmondGraphMatching.hpp"

using namespace boost;
using namespace Jaal;

typedef adjacency_list<vecS, vecS, undirectedS> Graph;

////////////////////////////////////////////////////////////////////////

int AllQuadMeshGenerator :: getEdmondsGraphMatching()
{
    if( trimesh == nullptr ) {
        cout << "Warning: Mesh object is null " << endl;
        return 1;
    }

    int topDim = trimesh->getTopology()->getDimension();
    if( topDim != 2 ) {
        cout << "Warning: Graph matching only for the triangle mesh " << endl;
        return 1;
    }

    size_t numFaces = trimesh->getActiveSize(2);

    if( numFaces%2) {
        cout << "Info: Creating even triangulation" << endl;
        AllQuadMeshGenerator::makeEvenTriangulation();
        numFaces = trimesh->getActiveSize(2);
        trimesh->pruneAll();
    }

    numFaces = trimesh->getActiveSize(2);
    assert( numFaces%2 == 0);

    trimesh->getTopology()->searchBoundary();

    JStopWatch swatch;

    swatch.start();

    JMeshDualGraph dgrapher;
    dgrapher.setMesh(trimesh);
    JMeshPtr dgraph = dgrapher.getGraph();
    int numnodes = dgraph->getSize(0);
    int numedges = dgraph->getSize(1);

    Graph graph(numnodes);
    for( int i = 0; i < numedges; i++) {
        const JEdgePtr &edge = dgraph->getEdgeAt(i);
        if( edge->isActive() ) {
            int v0 = edge->getNodeAt(0)->getID();
            int v1 = edge->getNodeAt(1)->getID();
            add_edge(v0,v1,graph);
        }
    }

    dgraph.reset();

    std::vector<graph_traits<Graph>::vertex_descriptor> mate(numnodes);

    bool success = checked_edmonds_maximum_cardinality_matching(graph, &mate[0]);
    assert(success);

    vector<JFacePair>  boostpairs;
    JFacePair newpair;
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for(tie(vi,vi_end) = vertices(graph); vi != vi_end; ++vi) {
        if (mate[*vi] != graph_traits<Graph>::null_vertex() && *vi < mate[*vi]) {
            newpair.first  = trimesh->getFaceAt( *vi );
            newpair.second = trimesh->getFaceAt( mate[*vi] );
            boostpairs.push_back( newpair );
        }
    }

    quadmesh = AllQuadMeshGenerator::collapse_matched_triangles(trimesh, boostpairs);
    swatch.stop();
    quadmesh->getTopology()->collectEdges();

    if( quadmesh->getTopology()->isHomogeneous(2) != JFace::QUADRILATERAL ) {
        cout << "Warning: All elements in the quadmesh are not quad " << endl;
    }

    if( !quadmesh->getTopology()->isConsistent() )
        quadmesh->getTopology()->getConsistent();

    if( quadmesh->getSize(2) != (size_t) (numFaces/2) )
        cout << "Warning: Imperfect matching has occured " << endl;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
JMeshPtr JEdmondGraphMatching :: getEdgeMatching()
{
    JMeshPtr wiremesh;

    if( mesh == nullptr) return wiremesh;

    wiremesh = JMesh::newObject();

    size_t numnodes = mesh->getSize(0);
    size_t numedges = mesh->getSize(1);

    JNodeSequence nodes = mesh->getNodes();
    wiremesh->addObjects(nodes);
    
    Graph graph(numnodes);
    for( int i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( !edge->isBoundary() ) {
                int v0 = edge->getNodeAt(0)->getID();
                int v1 = edge->getNodeAt(1)->getID();
                add_edge(v0,v1,graph);
            }
        }
    }

    std::vector<graph_traits<Graph>::vertex_descriptor> mate(numnodes);

    bool success = checked_edmonds_maximum_cardinality_matching(graph, &mate[0]);
    assert(success);

    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for(tie(vi,vi_end) = vertices(graph); vi != vi_end; ++vi) {
        if (mate[*vi] != graph_traits<Graph>::null_vertex() && *vi < mate[*vi]) {
            JNodePtr  v0 = mesh->getNodeAt(*vi);
            JNodePtr  v1 = mesh->getNodeAt(mate[*vi]);
            JEdgePtr  edge = JSimplex::getEdgeOf(v0, v1);
            wiremesh->addObject(edge);
        }
    }

    return wiremesh;
}
