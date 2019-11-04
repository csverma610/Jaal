#pragma once

#include "Curve.hpp"

struct Knot {
    std::vector<JCurve*> loops;
    std::vector<int>     getGaussCode();
    std::vector<int>     getDTCode();
    std::vector<int>     getConway();
};

class ClassicalKnots
{
public:
    Knot*  getKnot( int component, int crossing, int index );
    Knot*  getKnot( int id );
private:
    void readKnotFile();
};

