#include "MeshGeodesics.hpp"

void JMeshGeodesics :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    initialized = 0;
}

////////////////////////////////////////////////////////////////////////////////
void JTriMeshGeodesics :: initGoogleTriMesh( const JMeshPtr &triMesh)
{
    assert( triMesh != nullptr);
    // If we are using Google Geodescic software, all the elements must be
    // triangle ...
    size_t numnodes = triMesh->getActiveSize(0);
    points.resize(3*numnodes);
    size_t index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = triMesh->getNodeAt(i);
        if( vtx->isActive() ) {
            const Point3D &xyz = vtx->getXYZCoords();
            points[3*index+0] = xyz[0];
            points[3*index+1] = xyz[1];
            points[3*index+2] = xyz[2];
            index++;
        }
    }

    size_t numfaces = triMesh->getActiveSize(2);
    faces.resize(3*numfaces);
    index = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = triMesh->getFaceAt(i);
        if( face->isActive() ) {
            for( int j = 0; j < 3; j++)
                faces[3*index+j] = face->getNodeAt(j)->getID();
        }
    }
    geodesicTriMesh.reset( new geodesic::Mesh );
    geodesicTriMesh->initialize_mesh_data(points, faces);
}
////////////////////////////////////////////////////////////////////////////////

void JTriMeshGeodesics :: initMesh()
{
    if( mesh == nullptr) return;

    int topDim = mesh->getTopology()->getDimension();
    AllTriMeshGenerator alltri;
    if( topDim == 2 ) {
        int elemType = mesh->getTopology()->isHomogeneous(2);
        if( elemType == JFace::TRIANGLE)
            initGoogleTriMesh(mesh);

        if( algorithm == EXACT_DIJKSTRA && elemType == JFace::QUADRILATERAL) {
            JMeshPtr tmpmesh1 = mesh->deepCopy();
            JMeshPtr tmpmesh2 = alltri.fromQuadMesh(tmpmesh1, 4);
            initGoogleTriMesh(tmpmesh2);
        }
    }

    if( geodesicAlgorithm == nullptr) {
        if( algorithm == APPROXIMATE_DIJKSTRA)
            geodesicAlgorithm.reset( new geodesic::GeodesicAlgorithmDijkstra(geodesicTriMesh.get()));
        if( algorithm == EXACT_DIJKSTRA)
            geodesicAlgorithm.reset( new geodesic::GeodesicAlgorithmExact(geodesicTriMesh.get()));
        if( algorithm == SUBDIVISION_DIJKSTRA)
            geodesicAlgorithm.reset( new geodesic::GeodesicAlgorithmSubdivision(geodesicTriMesh.get()));
    }
    initialized = 1;
}

////////////////////////////////////////////////////////////////////////////////

vector<Point3D> JTriMeshGeodesics :: getExactPath( const JNodePtr &src, const JNodePtr &dst)
{
    vector<Point3D> pathCoords;
    geodesic::SurfacePoint source( &geodesicTriMesh->vertices()[src->getID()] );
    std::vector<geodesic::SurfacePoint> all_sources(1,source);
    geodesic::SurfacePoint target(&geodesicTriMesh->vertices()[dst->getID()]);

    std::vector<geodesic::SurfacePoint> path;
    geodesicAlgorithm->geodesic(source,target, path);
    if( !path.empty() ) {
        pathCoords.resize( path.size() );
        for(size_t i = 0; i<path.size(); ++i)
        {
            geodesic::SurfacePoint& s = path[i];
            pathCoords[i][0] =  s.x();
            pathCoords[i][1] =  s.y();
            pathCoords[i][2] =  s.z();
        }
    }
    return pathCoords;

}
///////////////////////////////////////////////////////////////////////////////////////////

JEdgeSequence JTriMeshGeodesics :: getApproxPath( const JNodePtr &src, const JNodePtr &dst)
{
    JEdgeSequence edges;
    vector<Point3D> pathCoords;
// getPath(src, dst, pathCoords);
    return edges;
}
///////////////////////////////////////////////////////////////////////////////////
int JTriMeshGeodesics :: setDistanceField( const JNodeSequence &srcnodes)
{
    if( !initialized) initMesh();

    geodesic::GeodesicAlgorithmExact  exactAlgo(geodesicTriMesh.get());
    std::vector<geodesic::SurfacePoint> sources( srcnodes.size() );

    size_t index = 0;
    for( const JNodePtr &vtx : srcnodes) {
        size_t id = vtx->getID();
        sources[index++] = geodesic::SurfacePoint(&geodesicTriMesh->vertices()[id]);
    }

    exactAlgo.propagate(sources);

    size_t numnodes = geodesicTriMesh->vertices().size();
    vector<double> dist( numnodes );
    double distance;
    for(size_t i=0; i< numnodes; ++i)
    {
        geodesic::SurfacePoint p(&(geodesicTriMesh->vertices())[i]);
        exactAlgo.best_source(p, distance);
        dist[i] = distance;
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////

void JGraphGeodesics :: initMesh()
{
    if( mesh == nullptr) return;

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &e = mesh->getEdgeAt(i);
        if( e->isActive() ) {
            double d =  JEdgeGeometry::getLength(e);
            e->setAttribute("Length", d);
        }
    }
    mesh->buildRelations(0,1);
    initialized = 1;
}

////////////////////////////////////////////////////////////////////////////////

void JGraphGeodesics :: initialize()
{
    double inftyDist = std::numeric_limits<double>::max();

    DistInfo  distInfo;
    distInfo.distance = 0.99*inftyDist;
    distInfo.prevNode = nullptr;

    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            distInfo.thisNode = v;
            v->setAttribute("DistInfo", distInfo);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

void JGraphGeodesics :: clear()
{
    if( mesh == nullptr) return;
    mesh->deleteNodeAttribute("DistInfo");
}

///////////////////////////////////////////////////////////////////////////////////////////
void JGraphGeodesics :: fastmarching( PriorityQ  &nodesQ, const JNodePtr &dst)
{
    double v0dist, v1dist, elen = 0.0;
    JEdgeSequence edges;
    DistInfo distInfo1,distInfo2;
    int err;

    while(!nodesQ.empty() ) {
        distInfo1 = nodesQ.top();
        JNodePtr currNode = distInfo1.thisNode;
        nodesQ.pop();
        JNode::getRelations(currNode, edges);
        assert( !edges.empty() ) ;
        v0dist = distInfo1.distance;
        for( size_t i = 0; i < edges.size(); i++) {
            const JNodePtr nextNode =  edges[i]->getOtherNode(currNode);
            err = edges[i]->getAttribute("Length", elen);
            assert(!err);
            err = nextNode->getAttribute("DistInfo", distInfo2);
            assert(!err);
            v1dist = distInfo2.distance;
            if( v1dist > v0dist + elen ) {
                distInfo2.distance = v0dist + elen;
                distInfo2.prevNode = currNode;
                nextNode->setAttribute("DistInfo", distInfo2);
                nodesQ.push(distInfo2);
            }
        }
        if( dst == nullptr) {
            distInfo2 = nodesQ.top();
            if( meshFilterPtr != nullptr) {
                if( !meshFilterPtr->passThrough(distInfo2.thisNode) ) return;
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////

int JGraphGeodesics :: traceback(const JNodePtr &src, const JNodePtr &dst, JNodeSequence &nodeSeq)
{
    nodeSeq.clear();

    if( src == nullptr ) return 1;
    if( dst == nullptr ) return 1;

    nodeSeq.push_back(dst);
    JNodePtr prevNode;
    DistInfo distInfo;

    JNodePtr currNode = dst;
    while(1) {
        currNode->getAttribute("DistInfo", distInfo);
        prevNode = distInfo.prevNode;
        if( prevNode == nullptr) break;
        nodeSeq.push_back(prevNode);
        currNode = prevNode;
    }
    std::reverse( nodeSeq.begin(), nodeSeq.end() );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////

JEdgeSequence JGraphGeodesics :: getPath(const JNodePtr &src, const JNodePtr &dst)
{
    JEdgeSequence edgeSeq;
    if( mesh == nullptr || src == nullptr ) return edgeSeq;

    if( initialized == 0) initMesh();

    PriorityQ  nodesQ;
    nodesQ = PriorityQ();

    DistInfo distInfo;
    distInfo.distance = 0.0;
    distInfo.thisNode = src;
    src->setAttribute("DistInfo", distInfo);
    nodesQ.push(distInfo);

    fastmarching(nodesQ , dst);

    JNodePtr  node2;
    if( dst == nullptr) {
        if( !nodesQ.empty() )
            distInfo = nodesQ.top();
        node2 = distInfo.thisNode;
    } else
        node2 = dst;

    JNodeSequence nodeSeq;
    traceback(src, node2, nodeSeq);

    if( nodeSeq.size()  < 2) return edgeSeq;

    size_t numEdges = nodeSeq.size()-1;
    edgeSeq.resize( numEdges );
    for( size_t i = 0; i < numEdges; i++) {
        const JNodePtr &v0  = nodeSeq[i];
        const JNodePtr &v1  = nodeSeq[i+1];
        JEdgePtr edge  = JSimplex::getEdgeOf(v0,v1);
        edgeSeq[i] = edge;
    }
    return edgeSeq;
}

///////////////////////////////////////////////////////////////////////////////////////////
int  JGraphGeodesics :: setDistanceField(const JNodeSequence &src)
{
    if( mesh == nullptr ) return 0;

    if( initialized == 0) initMesh();
    initialize();

    PriorityQ  nodesQ;
    nodesQ = PriorityQ();

    DistInfo distInfo;
    for( size_t i = 0; i < src.size(); i++) {
        src[i]->getAttribute("DistInfo", distInfo);
        distInfo.distance = 0.0;
        distInfo.thisNode = src[i];
        src[i]->setAttribute("DistInfo", distInfo);
        nodesQ.push(distInfo);
    }
    fastmarching(nodesQ , nullptr);

    assert( nodesQ.empty() );
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("DistInfo", distInfo);
            int id    = vtx->getID();
            vtx->setAttribute("Distance", distInfo.distance);
            vtx->deleteAttribute("DistInfo");
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
double JGraphGeodesics :: getDistance(const JNodePtr &src, const JNodePtr &dst)
{
    /*
        JEdgeSequence path;
        getApproxPath(src, dst, path);
        double len = JEdgeGeometry::getLength(path);
        return len;
    */
}
///////////////////////////////////////////////////////////////////////////////////////////
