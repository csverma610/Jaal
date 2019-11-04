#pragma once

#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

class JMeshHolesFill
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;
    typedef Polyhedron::Halfedge_handle    Halfedge_handle;
    typedef Polyhedron::Facet_handle       Facet_handle;
    typedef Polyhedron::Vertex_handle      Vertex_handle;

public:
    void setMesh( const JMeshPtr &m);
    int  getNumOfHoles() const;
    JMeshPtr fillAll();
private:
    JMeshPtr mesh;
    Polyhedron poly;
};
