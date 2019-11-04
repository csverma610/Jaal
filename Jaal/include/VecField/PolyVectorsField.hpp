#pragma once

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/n_polyvector.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "SurfaceVectorField.hpp"

class JPolyVectorsField : public JSurfaceVectorField
{
public:
    void setMesh( const JMeshPtr &m);
    void genRandomConstraints(int nRandom);
    int  genField();
};
