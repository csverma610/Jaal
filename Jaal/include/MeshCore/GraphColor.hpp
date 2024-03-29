#pragma once

#include "Mesh.hpp"

using namespace Jaal;

class JGraphColor
{
public:
    int  setColors(const JMeshPtr &m);
    int  getMaxColors() const
    {
        return maxIndex;
    }
private:
    JMeshPtr mesh;
    string name;
    int  maxIndex;
    JNodeSequence vneighs;
    vector<int> usedColors, unusedColors;

    void setColor( const JNodePtr &v);
    int  verify();
};

