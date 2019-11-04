#include "AllQuadMeshGenerator.hpp"


JMeshPtr hamiltonian_triangulation( const JMeshPtr &orgtrimesh )
{
    JMeshPtr newtrimesh = JMesh::newObject();

    size_t numnodes = orgtrimesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = orgtrimesh->getNodeAt(i);
        newtrimesh->addObject(vtx);
    }

    JFacePtr face, newtri;
    JNodePtr dualvtx;
    Point3D p3d;

    size_t numfaces = orgtrimesh->getSize(2);
    for( size_t iface = 0; iface < numfaces; iface++) {
        face = orgtrimesh->getFaceAt(iface);
        assert( face );
        face->setID(iface);
        face->getAvgXYZ( p3d );
        dualvtx = JNode::newObject();
        dualvtx->setID( numnodes + iface );
        dualvtx->setXYZCoords(p3d);
        face->setAttribute("Dual", dualvtx);
        newtrimesh->addObject(dualvtx);
    }

    orgtrimesh->buildRelations(1,2);

    JFaceSequence efaces;
    JNodePtr dual0, dual1;

    for( size_t iface = 0; iface < numfaces; iface++) {
        face = orgtrimesh->getFaceAt(iface);
        for( int j = 0; j < 3; j++) {
            JEdgePtr edge = face->getEdgeAt(j);
            JNodePtr v0 = edge->getNodeAt(0);
            JNodePtr v1 = edge->getNodeAt(1);
            JEdge::getRelations(edge, efaces);
            if( efaces.size() == 2 ) {
                if( min( efaces[0], efaces[1] ) == face ) {
                    efaces[0]->getAttribute("Dual", dual0 );
                    efaces[1]->getAttribute("Dual", dual1 );
                    JEdgePtr dedge = JEdge::newObject( dual0, dual1);
                    dedge->setAttribute("Dual", 1);
                    newtrimesh->addObject(edge);

                    newtri = Triangle::newObject( dual0, dual1, v0);
                    newtrimesh->addObject(newtri);

                    newtri = Triangle::newObject( dual0, dual1, v1);
                    newtrimesh->addObject(newtri);
                }
            } else {
                efaces[0]->getAttribute("Dual", dual0 );
                newtri = Triangle::newObject( dual0, v0, v1);
                newtrimesh->addObject(newtri);
            }

        }
    }

    orgtrimesh->deleteFaceAttribute("Dual");
    orgtrimesh->clearRelations(1,2);

    newtrimesh->getTopology()->getConsistent();

    return newtrimesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getHamiltonianQuads()
{
#ifdef CSV
    JMeshPtr trimesh = hamiltonian_triangulation( orgmesh );

    size_t numnodes = trimesh->getSize(0);

    JNodePtr v0, v1, ot1, ot2;

    JMeshPtr quadmesh = JMesh::newObject();

    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = trimesh->getNodeAt(i);
        quadmesh->addObject(vtx);
    }
//   trimesh->getTopology()->collect_edges();

    trimesh->buildRelations(1,2);

    JFaceSequence efaces;

    size_t numedges = trimesh->getSize(1);
    for( size_t iedge = 0; iedge < numedges; iedge++) {
        JEdgePtr edge = trimesh->getEdgeAt(iedge);
        if( !edge->hasAttribute("Dual") ) {
            JEdge::getRelations( edge, efaces );
            if( efaces.size() == 2 ) {
                if( efaces[0]->isActive() && efaces[1]->isActive() ) {
                    v0  = edge->getNodeAt(0);
                    v1  = edge->getNodeAt(1);
                    ot1 = Triangle::getOppositeNode( efaces[0], v0, v1);
                    ot2 = Triangle::getOppositeNode( efaces[1], v0, v1);
                    JFacePtr newquad =  Quadrilateral::newObject(ot1, v0, ot2, v1);
                    quadmesh->addObject(newquad);
                    efaces[0]->setStatus( JMeshEntity::REMOVE );
                    efaces[1]->setStatus( JMeshEntity::REMOVE );
                }
            }
        }
    }

    size_t numfaces = trimesh->getSize(2);
    for( size_t iface = 0; iface < numfaces; iface++) {
        JFacePtr face = trimesh->getFaceAt(iface);
        if( face->isActive() ) quadmesh->addObject( face );
    }

    trimesh->deleteEdgeAttribute("Dual");

    quadmesh->getTopology()->getConsistent();

    return quadmesh;

#ifdef CSV
    visitmark = v0->isVisited() + v1->isVisited();
    if( visitmark != 2) {
        neighs = Mesh::getRelation112( v0, v1 );
        if( neighs.size() == 2 ) {
            if( !neighs[0]->isRemoved() && !neighs[1]->isRemoved() ) {
                ot1 = Face::opposite_node( neighs[0], v0, v1);
                ot2 = Face::opposite_node( neighs[1], v0, v1);
                newquad =  new Face;
                connect[0] = ot1;
                connect[1] = v0;
                connect[2] = ot2;
                connect[3] = v1;
                newquad->setConnection(connect);
                quadmesh->addFace(newquad);
                neighs[0]->setRemoveMark(1);
                neighs[1]->setRemoveMark(1);
            }
        }
    }

    visitmark = v1->isVisited() + v2->isVisited();
    if( visitmark != 2) {
        neighs = Mesh::getRelation112( v1, v2 );
        if( neighs.size() == 2 ) {
            if( !neighs[0]->isRemoved() && !neighs[1]->isRemoved() ) {
                ot1 = Face::opposite_node( neighs[0], v1, v2);
                ot2 = Face::opposite_node( neighs[1], v1, v2);
                newquad =  new Face;
                connect[0] = ot1;
                connect[1] = v1;
                connect[2] = ot2;
                connect[3] = v2;
                newquad->setConnection(connect);
                quadmesh->addFace(newquad);
                neighs[0]->setRemoveMark(1);
                neighs[1]->setRemoveMark(1);
            }
        }
    }

    visitmark = v2->isVisited() + v0->isVisited();
    if( visitmark != 2) {
        neighs = Mesh::getRelation112( v2, v0 );
        if( neighs.size() == 2 ) {
            if( !neighs[0]->isRemoved() && !neighs[1]->isRemoved() ) {
                ot1 = Face::opposite_node( neighs[0], v2, v0);
                ot2 = Face::opposite_node( neighs[1], v2, v0);
                newquad =  new Face;
                connect[0] = ot1;
                connect[1] = v2;
                connect[2] = ot2;
                connect[3] = v0;
                newquad->setConnection(connect);
                quadmesh->addFace(newquad);
                neighs[0]->setRemoveMark(1);
                neighs[1]->setRemoveMark(1);
            }
        }
    }
#endif
}

for( int iface = 0; iface < numfaces; iface++) {
    face = trimesh->getFace(iface);
    if( !face->isRemoved() ) quadmesh->addFace(face);
}

if( !relexist) trimesh->clearRelations(0,2);

numfaces = quadmesh->getSize(2);
int numtris = 0, numquads = 0;
for( int i = 0; i < numfaces; i++) {
    int nnodes = quadmesh->getFace(i)->getSize(0);
    if( nnodes == 3 ) numtris++;
    if( nnodes == 4 ) numquads++;
}
quadmesh->enumerate(2);

cout << "# Triangles : " << numtris  << endl;
cout << "# Quads     : " << numquads << endl;

QuadCleanUp quadclean(quadmesh);

quadclean.search_diamonds();
quadclean.search_doublets();

return quadmesh;
#endif
return nullptr;

}

///////////////////////////////////////////////////////////////////////////////
