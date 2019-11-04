#pragma once
#include "Mesh.hpp"

struct TriSimplify {
    int simplify( Face *f, JNodeSequence &newnodes, JFaceSequence &newfaces);
    int simplify( Mesh *m, size_t numfaces);
};
