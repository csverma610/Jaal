#include "MeshPartitioner.hpp"

//////////////////////////////////////////////////////////////////////////////////

int JMetisPartitioner :: getPartitions( int n)
{
    if( mesh == nullptr) return 1;
    numParts = n;

    clear();

    mesh->enumerate(0);
    int topDim = mesh->getTopology()->getDimension();

    vector<size_t> eConnect;
    vector<idx_t> eind, eptr, epart, npart;

    idx_t nn = mesh->getSize(0);
    npart.resize(nn);

    int elmtype = mesh->getTopology()->isHomogeneous(topDim);

    vector<int> topo; // Not used ...
    mesh->getTopology()->getNodesArray( eConnect, topo);

    eind.resize(eConnect.size());

    boost::copy( eConnect, eind.begin() );

    idx_t ne = mesh->getSize(topDim);
    epart.resize(ne);
    eptr.resize(ne+1);

    int index;
    int etype = 0;
    if( topDim == 2 ) {
        switch(elmtype) {
        case JFace::TRIANGLE:
            etype = 1;
            break;
        case JFace::QUADRILATERAL:
            etype = 4;
            break;
        }
        eptr[0]  = 0;
        index    = 1;
        int nsum = 0;
        for( idx_t i = 0; i < ne; i++)  {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                nsum += f->getSize(0);
                eptr[index++] = nsum;
            }
        }
    }

    if( topDim == 3 ) {
        switch(elmtype) {
        case JCell::TETRAHEDRON:
            etype = 2;
            break;
        case JCell::HEXAHEDRON:
            etype = 3;
            break;
        }
        eptr[0]  = 0;
        index    = 1;
        int nsum = 0;
        for( idx_t i = 0; i < ne; i++)  {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() ) {
                nsum += c->getSize(0);
                eptr[index++] = nsum;
            }
        }
    }

    if( etype  == 0) return 1;

#ifdef METIS_VERSION_LESS_THAN_5
    int    numflag = 0, edgecut;
    METIS_PartMeshDual(&ne, &nn, &eind[0], &etype, &numflag, &numParts,
                       &edgecut, &epart[0], &npart[0]);
#else
    idx_t ncommon = 2;
    idx_t objval;
    idx_t nparts = numParts;
    METIS_PartMeshDual(&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval,
                       &epart[0], &npart[0]);
#endif

    if( topDim == 2 ) {
        index = 0;
        for( int i = 0; i < ne; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) face->setAttribute("Partition", (int)epart[index++] );
        }
    }

    if( topDim == 3 ) {
        index = 0;
        for( int  i = 0; i < ne; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) cell->setAttribute("Partition", (int)epart[index++] );
        }
    }

    index = 0;
    for( int i = 0; i < nn; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) vertex->setAttribute("Partition", (int)npart[index++] );
    }

    searchInterfaces();
    searchCorners();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
void JMetisPartitioner2D::split_disjoint_component(int cid)
{
    /*
        SharedFaceMesh partmesh = getPart(cid);

        if (partmesh->getSize() == 0) return;

        SharedFace face;
        deque<SharedVertex> nodes_to_visit;

        face = partmesh->at(0);
        for (int i = 0; i < face->getNumNodes(); i++)
            nodes_to_visit.push_back(face->getVertex(i));

        partmesh->buildRelations(0, 2);

        FOR_EACH_FACE(face, partmesh)
             face->setVisit(0);

        FaceContainer neighs;
        SharedVertex curr_vertex;

        while (!nodes_to_visit.empty())
        {
            curr_vertex = nodes_to_visit.front();
            nodes_to_visit.pop_front();
            neighs = curr_vertex->getRelations02();
            for (int i = 0; i < neighs.size(); i++)
            {
                face = neighs[i];
                if (!face->getVisit())
                {
                    for (int j = 0; j < face->getNumNodes(); j++)
                        nodes_to_visit.push_back(face->getVertex(j));
                }
                face->setVisit(1);
            }
        }

        int index = 0;
        FOR_EACH_FACE(face, partmesh)
        {
            if (!face->getVisit())
            {
                face->setAttribute(eKey, numParts);
                index++;
            }
        }
        if (index) numParts++;

        partmesh->clearRelations(0, 2);
    */
}

///////////////////////////////////////////////////////////////////////////////

#endif

