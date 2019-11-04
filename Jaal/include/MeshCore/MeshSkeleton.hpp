#pragma once

#include "Mesh.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
//#include <CGAL/Mean_curvature_flow_skeletonization.h>


using namespace std;

class JMeshSkeleton
{
    typedef CGAL::Simple_cartesian<double>                        Kernel;
    typedef Kernel::Point_3                                       Point;
    typedef CGAL::Surface_mesh<Point>                             Triangle_mesh;
    /*
        typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
        typedef Skeletonization::Skeleton                             Skeleton;
        typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
        typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
        typedef Skeleton::edge_descriptor                             Skeleton_edge;
    */
public:
    void setMesh(const JMeshPtr &m)
    {
        mesh = m;
    }

    JNodeSequence getMedialNodes() const;

private:
    JMeshPtr mesh;
};

