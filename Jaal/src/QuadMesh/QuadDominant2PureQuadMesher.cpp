#include "QuadDominant2PureQuadsMesher.hpp"

void JQuadDominant2PureQuadsMesher :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    firstNonQuadFilter.reset( new FirstNonQuad );
}

/////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsMesher :: setOnlyTrianglesAsNonQuads()
{
    if( mesh == nullptr ) return;

    size_t numfaces = mesh->getSize(2);

    int err;
    JNodeSequence  connect;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive()) {
            int nn = face->getSize(0);
            if( nn == 5) {
                connect = face->getNodes();
                JFacePtr quad = Quadrilateral::newObject( connect[0], connect[1], connect[2], connect[3] );
                JFacePtr tri = Triangle::newObject( connect[0], connect[3], connect[4]);
                err = mesh->addObject(quad);
                assert(!err);
                err = mesh->addObject(tri);
                assert(!err);
                face->setStatus(JMeshEntity::REMOVE);
            }
        }
    }
    mesh->pruneFaces();
}
/////////////////////////////////////////////////////////////////////////////
size_t JQuadDominant2PureQuadsMesher :: getNumOfNonQuads() const
{
    if( mesh == nullptr ) return 0;

    size_t numfaces = mesh->getSize(2);

    size_t nCount = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) != 4) nCount++;
    }

    return nCount;
}

/////////////////////////////////////////////////////////////////////////////
vector<JFaceSequence> JQuadDominant2PureQuadsMesher :: getAllStrips()
{
    quadStrips.clear();

    vector<JFaceSequence> emptyQ;
    if( mesh == nullptr ) return emptyQ;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() )
            face->setVisitBit(0);
    }

    set<JFacePtr> nonQuads;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) != 4)
            nonQuads.insert(face);
    }

    if( nonQuads.empty() ) return emptyQ;

    JMeshDualGraph  dg;
    dg.setMesh(mesh);
    JMeshPtr dualGraph = dg.getGraph();

    JGraphGeodesics  geode;
    geode.setMesh( dualGraph );
    geode.setFilter(firstNonQuadFilter);

    size_t numStrips = nonQuads.size()/2;

    JFacePtr startFace, endFace, pface;
    JNodePtr dnode;
    JFaceSequence newStrip;
    JNodeSequence dualnodes;
    for( size_t i = 0; i < numStrips; i++) {
        geode.initialize();
        startFace = *nonQuads.begin();
        nonQuads.erase(startFace);
        startFace->getAttribute("DualNode", dnode);
        JEdgeSequence path = geode.getPath(dnode);
        JEdgeTopology::getChain(path, dnode);
        JNodePtr endNode = path.back()->getNodeAt(1);
        endNode->getAttribute("PrimalFace", endFace);
        if( nonQuads.find(endFace) != nonQuads.end() ) {
            nonQuads.erase(endFace);
            JEdgeTopology::getChainNodes(path, dualnodes);
            newStrip.clear();
            for( const JNodePtr &vtx: dualnodes) {
                vtx->getAttribute("PrimalFace", pface);
                newStrip.push_back(pface);
            }
            startFace->setVisitBit(1);
            endFace->setVisitBit(1);
            quadStrips.push_back(newStrip);
        }
    }

    boost::sort(quadStrips, []( const JFaceSequence &a, const JFaceSequence &b)
    {
        return a.size() < b.size();
    });

    numStrips = quadStrips.size();
    vector<JFaceSequence> vq(numStrips);
    for( int i = 0; i < numStrips; i++)
        vq[i] = quadStrips[i];

    return vq;
}

/////////////////////////////////////////////////////////////////////////////
JFaceSequence JQuadDominant2PureQuadsMesher :: getCurrentStrip()
{
    if( quadStrips.empty() ) getAllStrips();

    JFaceSequence emptyStrip;
    if( quadStrips.empty() ) return emptyStrip;

    size_t nstrips = quadStrips.size();
    for( size_t i = 0; i < nstrips; i++) {
        if( quadStrips.empty() ) emptyStrip;
        bool active = 1;
        for( const JFacePtr &f: quadStrips[0]) {
            if( !f->isActive() )
                active = 0;
            break;
        }
        if( active) return quadStrips[0];
        quadStrips.pop_front();
    }
    return emptyStrip;
}

//////////////////////////////////////////////////////////////////////////////
JFaceSequence JQuadDominant2PureQuadsMesher :: getRefineStrip()
{
}

//////////////////////////////////////////////////////////////////////////////

int JQuadDominant2PureQuadsMesher :: flipCurrentStrip()
{
    if( quadStrips.empty() ) return 1;

    JFaceSequence currSeq = quadStrips[0];
    quadStrips.pop_front();

    int nsize = currSeq.size();
    if( nsize < 2) {
        cout << "Warning: At least two faces are required in the sequence " << endl;
        return 1;
    }

    if( currSeq.front()->getSize(0) != 3 ) {
        cout << "Warning: Begining of the sequence is not a triangle " << endl;
        return 1;
    }

    if( currSeq.back()->getSize(0) != 3 ) {
        cout << "Warning: End of the sequence is not a triangle " << endl;
        return 1;
    }

    for( int i = 1; i < nsize-1; i++) {
        if( currSeq[i]->getSize(0) != 4) {
            cout << "Warning: Intermediate elements must be quad " << endl;
            return 1;
        }
    }

    JEdgePtr commonEdge, nextEdge, newEdge;
    JEdgeSequence commonEdges;
    JNodeSequence nodes(5);

    deque<JFacePtr> faceQ;
    for( const JFacePtr &f : currSeq) faceQ.push_back(f);

    JFacePtr newQuad, newTri;

    for( int i = 0; i < nsize-2; i++) {
        const JFacePtr &tri = faceQ[0];        // Must be triangle ...
        const JFacePtr &quad = faceQ[1];       // Must be Quad ..
        const JFacePtr &nextFace = faceQ[2];   // Can be any face ...
        faceQ.pop_front();
        faceQ.pop_front();

        // First face must be triangle and second facemust be quad.
        if( (tri->getSize(0) != 3) && (quad->getSize(0) != 4) ) break;
        JFace::getSharedEntities( tri, quad, commonEdges);
        if( commonEdges.size() != 1 ) return 1;
        commonEdge = commonEdges[0];

        // With respect to the quad, what is the position of "commonEdge"
        int pos0 = quad->getPosOf(commonEdge);
        assert( pos0 >= 0 && pos0 < 4);

        // Which edge the Quad and Next face sharing ...
        JFace::getSharedEntities( quad, nextFace, commonEdges);
        if( commonEdges.size() != 1 ) return 1;
        nextEdge = commonEdges[0];

        // What is the position of "NextEdge" relative to pos0 ? It can be 1,2, or 3
        int pos1 = -1;
        if( quad->getEdgeAt(pos0+1) == nextEdge ) pos1 = 1;
        if( quad->getEdgeAt(pos0+2) == nextEdge ) pos1 = 2;
        if( quad->getEdgeAt(pos0+3) == nextEdge ) pos1 = 3;

        assert( pos1 >= 0);

        nodes[0] = commonEdge->getNodeAt(0);
        nodes[2] = commonEdge->getNodeAt(1);
        nodes[1] = Triangle::getOppositeNode( tri, nodes[0], nodes[2] );
        int pos3 = tri->getPosOf( nodes[1] );
        nodes[2] = tri->getNodeAt(pos3+1);
        nodes[0] = tri->getNodeAt(pos3+2);
        nodes[3] = Quadrilateral::getDiagonalNode( quad, nodes[0] );
        nodes[4] = Quadrilateral::getDiagonalNode( quad, nodes[2] );

        switch(pos1)
        {
        case 1:
            newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[4]);
            newTri  = Triangle::newObject( nodes[2], nodes[3], nodes[4] );
            newEdge = JSimplex::getEdgeOf( nodes[2], nodes[4],1);
            break;
        case 2:
            newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[4]);
            newTri  = Triangle::newObject( nodes[2], nodes[3], nodes[4] );
            newEdge = JSimplex::getEdgeOf( nodes[2], nodes[4],1);
            break;
        case 3:
            newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[3]);
            newTri  = Triangle::newObject( nodes[0], nodes[3], nodes[4] );
            newEdge = JSimplex::getEdgeOf( nodes[0], nodes[3],1);
            break;
        default:
            return 1;
        }

        faceQ.push_front(newTri);

        mesh->addObject(newEdge);
        mesh->addObject(newTri);
        mesh->addObject(newQuad);

        tri->setStatus( JMeshEntity::REMOVE);
        quad->setStatus( JMeshEntity::REMOVE);
        commonEdge->setStatus( JMeshEntity::REMOVE);
    }

    assert( faceQ.size() == 2);
    assert( faceQ[0]->getSize(0) == 3);
    assert( faceQ[1]->getSize(0) == 3);

    JFace::getSharedEntities( faceQ[0], faceQ[1], commonEdges);
    assert( commonEdges.size() == 1);

    nodes[0] = commonEdges[0]->getNodeAt(0);
    nodes[2] = commonEdges[0]->getNodeAt(1);
    nodes[1] = Triangle::getOppositeNode( faceQ[0], nodes[0], nodes[2] );
    nodes[3] = Triangle::getOppositeNode( faceQ[1], nodes[0], nodes[2] );

    newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[3]);
    mesh->addObject(newQuad);
    faceQ[0]->setStatus( JMeshEntity::REMOVE);
    faceQ[1]->setStatus( JMeshEntity::REMOVE);
    commonEdges[0]->setStatus( JMeshEntity::REMOVE);

    mesh->pruneFaces();
    mesh->pruneEdges();
    return 0;
}
//////////////////////////////////////////////////////////////////////////////
