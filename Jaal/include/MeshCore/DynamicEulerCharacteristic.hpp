#pragma once

#include "Mesh.hpp"

class JDynamicEulerCharacteristic
{
public:
    void reset();
    int  addObject(const JMeshPtr &m);

    int addObject( const JCellPtr &);
    int addObject( const JFacePtr &);
    int addObject( const JEdgePtr &);
    int addObject( const JNodePtr &);

    int removeObject( const JCellPtr &);
    int removeObject( const JFacePtr &);
    int removeObject( const JEdgePtr &);
    int removeObject( const JNodePtr &);

    vector<size_t>  getFVector() const;

    int  getEulerCharacteristic() const
    {
        return nodeSet.size()-edgeSet.size()+ faceSet.size() - cellSet.size();
    }

private:
    long  birth_time;
    std::unordered_set<JNodePtr>  nodeSet;
    std::unordered_set<JEdgePtr>  edgeSet;
    std::unordered_set<JFacePtr>  faceSet;
    std::unordered_set<JCellPtr>  cellSet;
};

