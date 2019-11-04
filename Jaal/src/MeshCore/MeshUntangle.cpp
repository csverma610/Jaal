#include "MeshUntangle.hpp"
#include "StopWatch.hpp"

////////////////////////////////////////////////////////////////////////

JMeshUntangle :: JMeshUntangle()
{
    energyType = JLocallyInjectiveMap::DIRICHLET_ENERGY;
}

JMeshUntangle :: ~JMeshUntangle()
{
    if( inflatedMesh == nullptr ) return;
    triMesh.reset();
}

///////////////////////////////////////////////////////////////////////////////////
void JMeshUntangle :: setMesh( const JMeshPtr &m)
{
    srcMesh = m;
    if( srcMesh == nullptr) return;

    srcMesh->pruneAll();
    srcMesh->enumerate(0);

    srcMesh->getTopology()->getConsistent();

    entityDim = srcMesh->getTopology()->getDimension();

    // If all the elements are positive, we do not have to do anything.
    size_t nCount = countInverted( srcMesh );

//  if( nCount == 0) return;
    // Make sure that the mesh is topological consistent
    //  srcMesh->getTopology()->getConsistent();

    // If more than half are inverted, then flip the mesh ...
    if( entityDim == 2) {
        if( nCount > 0.5*srcMesh->getSize(2) ) {
            size_t numfaces = srcMesh->getSize(2);
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &f = srcMesh->getFaceAt(i);
                f->reverse();
            }
        }

        // For untangling purposes, we will be working with a triangle
        // mesh so make a trimesh from a quad mesh...
        if( srcMesh->getTopology()->isHomogeneous(2) == JFace::TRIANGLE) {
            triMesh = srcMesh->deepCopy();
        } else {
            JMeshPtr tmpmesh = srcMesh->deepCopy();
            triMesh = AllTriMeshGenerator::fromQuadMesh( tmpmesh, 4);
        }
    }

    if( entityDim == 3) {
        if( nCount > 0.5*srcMesh->getSize(3) ) {
            size_t numcells = srcMesh->getSize(3);
            for( size_t i = 0; i < numcells; i++) {
                const JCellPtr &c = srcMesh->getCellAt(i);
                c->reverse();
            }
        }

        // For untangling purposes, we will be working with a triangle
        // mesh so make a trimesh from a quad mesh...
        if( srcMesh->getTopology()->isHomogeneous(3) == JCell::HEXAHEDRON) {
            JMeshPtr tmpmesh = srcMesh->deepCopy();
            AllTetMeshGenerator  alltet;
            tetMesh = alltet.fromHexMesh(tmpmesh);
        } else
            tetMesh = srcMesh->deepCopy();
    }
}

////////////////////////////////////////////////////////////////////////
int JMeshUntangle :: smoothCorner(const JNodePtr &vtx)
{
    JNodeSequence vneighs;
    vtx->getRelations(vtx, vneighs);

    Point3D pmid;
    int nCount = 0;

    pmid[0] = 0.0;
    pmid[1] = 0.0;
    pmid[2] = 0.0;
    for( const JNodePtr &v : vneighs) {
        if( v->isBoundary() ) {
            const Point3D &p = v->getXYZCoords();
            pmid[0] += p[0];
            pmid[1] += p[1];
            pmid[2] += p[2];
            nCount++;
        }
    }
    if( nCount == 0) return 1;
    pmid[0] /= ( double)nCount;
    pmid[1] /= ( double)nCount;
    pmid[2] /= ( double)nCount;
    vtx->setXYZCoords(pmid);
    return 0;
}
////////////////////////////////////////////////////////////////////////
int JMeshUntangle :: smoothBoundary()
{
    if( inflateAll ) {
        size_t numnodes = triMesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = triMesh->getNodeAt(i);
            if( vtx->isBoundary() ) smoothCorner( vtx);
        }
        return 0;
    }

    for( const JNodePtr &vtx : concaveCorners)
        smoothCorner( vtx);

    return 0;
}

////////////////////////////////////////////////////////////////////////
JFaceSequence JMeshUntangle :: getInvertedFaces()
{
    JFaceSequence inverted;
    size_t numfaces = srcMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = srcMesh->getFaceAt(i);
        if( JFaceGeometry::isInverted(face) ) inverted.push_back(face);
    }
    return inverted;
}
////////////////////////////////////////////////////////////////////////
JCellSequence JMeshUntangle :: getInvertedCells()
{
    JCellSequence inverted;

    size_t numcells = srcMesh->getSize(3);
    for( size_t i = 0; i < numcells; i++) {
        const JCellPtr &cell = srcMesh->getCellAt(i);
        if( JCellGeometry::isInverted(cell) ) inverted.push_back(cell);
    }
    return inverted;
}
////////////////////////////////////////////////////////////////////////

size_t JMeshUntangle :: countInverted( const JMeshPtr &m)
{
    if( m == nullptr) return 0;
    concaveCorners.clear();

    size_t nCount = 0;

    if( entityDim == 2 ) {
        size_t numfaces = m->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = m->getFaceAt(i);
            if( JFaceGeometry::isInverted(face) ) {
                for( int j = 0; j < face->getSize(0); j++) {
                    const JNodePtr &v = face->getNodeAt(j);
                    if( v->isBoundary() ) {
                        concaveCorners.insert(v);
                    }
                }
                nCount++;
            }
        }
        return nCount;
    }

    if( entityDim == 3 ) {
        size_t numcells = m->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = m->getCellAt(i);
            if( JCellGeometry::isInverted(cell) ) {
                for( int j = 0; j < cell->getSize(0); j++) {
                    const JNodePtr &v = cell->getNodeAt(j);
                    if( v->isBoundary() ) concaveCorners.insert(v);
                }
                nCount++;
            }
        }
        return nCount;
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////
size_t JMeshUntangle :: getNumNodesConvexified() const
{
    if(inflatedMesh == nullptr) return 0;

}
////////////////////////////////////////////////////////////////////////

void JMeshUntangle :: backProject()
{
    JStopWatch swatch;
    JLocallyInjectiveMap lim;

    if( entityDim == 2) {
        if( triMesh == nullptr) return;
        triMesh->deleteNodeAttribute("Constraint");

        lim.setMesh(triMesh);
        lim.setEnergyType( energyType );
        lim.setMaxIterations( maxBackProjections);
        lim.solve();

        // Update the source ....
        size_t numnodes = srcMesh->getSize(0);
        #pragma omp parallel for
        for( size_t i = 0; i < numnodes; i++) {
            const Point3D &p = triMesh->getNodeAt(i)->getXYZCoords();
            const JNodePtr &vtx = srcMesh->getNodeAt(i);
            vtx->setXYZCoords(p);
        }
        JMeshGeometry jm(srcMesh);
        area[0] = jm.getArea();
    }

    if( entityDim == 3) {
        if( tetMesh == nullptr) return;
        tetMesh->deleteNodeAttribute("Constraint");

        lim.setMesh(tetMesh);
        lim.solve();

        // Update the source ....
        size_t numnodes = srcMesh->getSize(0);
        #pragma omp parallel for
        for( size_t i = 0; i < numnodes; i++) {
            const Point3D &p = tetMesh->getNodeAt(i)->getXYZCoords();
            const JNodePtr &vtx = srcMesh->getNodeAt(i);
            vtx->setXYZCoords(p);
        }
        tetMesh.reset();
    }

    if( countInverted(srcMesh) != 0)
        cout << "Sorry: The mesh could not be inverted " << endl;
}

////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshUntangle :: getInflatedMesh()
{
    if( srcMesh == nullptr) return nullptr;

    if( entityDim == 2) return getInflated2D();
    if( entityDim == 3) return getInflated3D();

    return inflatedMesh;
}

////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshUntangle :: getSimplicialMesh()
{
    if( srcMesh == nullptr) return nullptr;

    if( entityDim == 2) return triMesh;
    if( entityDim == 3) return tetMesh;

    return inflatedMesh;
}

////////////////////////////////////////////////////////////////////////
double JMeshUntangle :: getMaxDistance() const
{
    Point3D p0, p1;

    double maxDistance = 0.0;
    size_t numnodes = triMesh->getSize(0);
//#pragma omp parallel fohhhhhhhhhhhh
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = triMesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int err = vtx->getAttribute("TargetPos", p0);
            if( !err ) {
                const Point3D  &p1  = vtx->getXYZCoords();
                double d = JMath::length(p0, p1);
                maxDistance = max(maxDistance, d);
            }
        }
    }
    return maxDistance;
}
////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshUntangle :: getInflated2D()
{
    if( triMesh == nullptr) return nullptr;

    int nstart = 0;
    size_t numnodes;

    inflatedMesh = srcMesh->deepCopy();
    int niter = countInverted( inflatedMesh );
    if( niter == 0) return nullptr;

    triMesh->getTopology()->searchBoundary();
    numnodes = triMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = triMesh->getNodeAt(i);
        if( v->isBoundary() ) {
            const Point3D  &p  = v->getXYZCoords();
            v->setAttribute("TargetPos", p);
        }
    }

    triMesh->buildRelations(0,0);
    lap.setMesh(triMesh);
    lap.setNumIterations(100);
    lap.setBoundaryPreserve(1);
    lap.smoothAll();
    countInverted(triMesh);

//  niter = nstart;
    cout << "#Inverted Starts " << nstart << endl;
    int miter = 0;
    while(1) {
        smoothBoundary();
        lap.smoothAll();
        int nCount = countInverted( triMesh );
        cout << "#Inverted triangles "  << nCount << endl;
        if( nCount == 0 ||  miter++ >= maxBoundaryModifications) break;
    }

    // Copy the nodes from the triangle mesh to the "inflated mesh"
    numnodes = srcMesh->getSize(0);
    offset.resize(numnodes);
//#pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = inflatedMesh->getNodeAt(i);
        const JNodePtr &v1 = triMesh->getNodeAt(i);
        const Point3D  &p1 = v1->getXYZCoords();
        v0->setXYZCoords(p1);
        offset[i] = JNodeGeometry::getLength(srcMesh->getNodeAt(i), triMesh->getNodeAt(i));
    }
    cout << "#Inverted on input mesh    " << countInverted(srcMesh) << endl;
    cout << "#Inverted on triangle mesh " << countInverted(triMesh) << endl;
    cout << "#Inverted on inflated mesh " << countInverted(inflatedMesh) << endl;

    JMeshGeometry jm(inflatedMesh);
    area[1] = jm.getArea();

    return inflatedMesh;
}

////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshUntangle :: getInflated3D()
{
    if( tetMesh == nullptr) return nullptr;

    int nstart = 0;
    size_t numnodes;

    inflatedMesh = srcMesh->deepCopy();
    int niter = countInverted( inflatedMesh );
    if( niter == 0) return nullptr;

    tetMesh->getTopology()->searchBoundary();
    numnodes = tetMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = tetMesh->getNodeAt(i);
        if( v->isBoundary() ) {
            const Point3D  &p  = v->getXYZCoords();
            v->setAttribute("TargetPos", p);
        }
    }

    tetMesh->buildRelations(0,0);
    lap.setMesh(tetMesh);
    lap.setNumIterations(100);
    lap.setBoundaryPreserve(1);
    lap.smoothAll();
    countInverted(tetMesh);

//  niter = nstart;
    cout << "#Inverted Starts " << niter << endl;
    for( int i = 0; i < niter ; i++) {
        smoothBoundary();
        lap.smoothAll();
        int nCount = countInverted( tetMesh );
        cout << "#Inverted Tets: "  << nCount << endl;
        /*
                if( i%10 == 0) {
                    if( nCount >= nstart) {
                        cout << "Warning: Smoothing not reducing inverted elements " << endl;
                        break;
                    }
                    nstart = nCount;
                }
        */
        if( nCount == 0) break;
    }

    // Copy the nodes from the triangle mesh to the "inflated mesh"
    numnodes = srcMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = inflatedMesh->getNodeAt(i);
        const JNodePtr &v1 = tetMesh->getNodeAt(i);
        const Point3D  &p1 = v1->getXYZCoords();
        v0->setXYZCoords(p1);
    }
    cout << "#Inverted on input mesh    " << countInverted(srcMesh) << endl;
    cout << "#Inverted on tet   mesh    " << countInverted(tetMesh) << endl;
    cout << "#Inverted on inflated mesh " << countInverted(inflatedMesh) << endl;

    return inflatedMesh;
}

////////////////////////////////////////////////////////////////////////


