#pragma once

#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

class JMeshSDF
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> CGALPolyhedron;
public:
    void setMesh(const JMeshPtr &m)
    {
        mesh = m;
    }

    int  getNumOfSegments() const;
    int  execute();
private:
    JMeshPtr mesh;
};
