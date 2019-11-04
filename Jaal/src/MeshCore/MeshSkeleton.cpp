#include "MeshSkeleton.hpp"

JNodeSequence JMeshSkeleton :: getMedialNodes() const
{
    /*
       JNodeSequence midnodes;
       if( mesh == nullptr) return midnodes;

      JMeshIO meshio;
      meshio.saveAs(mesh, "tmp.off");
      std::ifstream input("tmp.off");

      Triangle_mesh tmesh;
      input >> tmesh;

      Skeleton skeleton;
      Skeletonization mcs(tmesh);
      // 1. Contract the mesh by mean curvature flow.
      mcs.contract_geometry();
      // 2. Collapse short edges and split bad triangles.
      mcs.collapse_edges();
      mcs.split_faces();
      // 3. Fix degenerate vertices.
      mcs.detect_degeneracies();
      // Perform the above three steps in one iteration.
      mcs.contract();
      // Iteratively apply step 1 to 3 until convergence.
      mcs.contract_until_convergence();
      // Convert the contracted mesh into a curve skeleton and
      // get the correspondent surface points
      mcs.convert_to_skeleton(skeleton);
      std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
      std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
      // Output all the edges of the skeleton.

      Point3D xyz;
      BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
         xyz[0] = skeleton[v].point.x();
         xyz[1] = skeleton[v].point.x();
         xyz[2] = skeleton[v].point.x();
      }
    */

    /*
      std::ofstream output("skel.cgal");
      BOOST_FOREACH(Skeleton_edge e, edges(skeleton))
      {
        const Point& s = skeleton[source(e, skeleton)].point;
        const Point& t = skeleton[target(e, skeleton)].point;
        output << "2 "<< s << " " << t << "\n";
      }
      output.close();

      // Output skeleton points and the corresponding surface points
      output.open("correspondance.cgal");
      BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
        BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
          output << "2 " << skeleton[v].point << "  " << get(CGAL::vertex_point, tmesh, vd)  << "\n";
    */
}






