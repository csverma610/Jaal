#pragma once

#include "Mesh.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "MeshGeodesics.hpp"

class JQuadDominant2PureQuadsMesher
{
    struct FirstNonQuad : public JMeshFilter {
        bool passThrough( const JNodePtr &v) const
        {
            if( v == nullptr) return 1;
            JFacePtr pface;
            int err = v->getAttribute("PrimalFace", pface);
            if( pface->getSize(0) == 4) return 1;
            if( pface->getVisitBit() == 1) return 1;
            return 0;
        }
    };

public:
    JQuadDominant2PureQuadsMesher() { }

    void setMesh(const JMeshPtr &m);
    void setOnlyTrianglesAsNonQuads();

    size_t getNumOfNonQuads() const;

    int mergeTriangles();
    int minRefinement();

    vector<JFaceSequence> getAllStrips();

    JFaceSequence getCurrentStrip();
    int  flipCurrentStrip();

    JFaceSequence getRefineStrip();
    int  refineCurrentStrip();

private:
    JMeshPtr mesh;
    std::shared_ptr<FirstNonQuad> firstNonQuadFilter;
    std::scoped_ptr<JMeshGeodesics> mGeodesics;
    deque<JFaceSequence> quadStrips;
    void  mergeTriQuad( const JFacePtr &tri, const JFacePtr &quad);
    void  genDualGraph();
};
