#include "MeshHolesFill.hpp"

JMeshPtr JMeshHolesFill::fillAll()
{
    size_t nb_holes = 0;
    BOOST_FOREACH(Halfedge_handle h, halfedges(poly))
    {
        if(h->is_border())
        {
            std::vector<Facet_handle>  patch_facets;
            std::vector<Vertex_handle> patch_vertices;
            bool success = CGAL::cpp11::get<0>(
                               CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                                   poly,
                                   h,
                                   std::back_inserter(patch_facets),
                                   std::back_inserter(patch_vertices),
                                   CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, poly)).
                                   geom_traits(Kernel())) );
            std::cout << " Number of facets in constructed patch: " << patch_facets.size() << std::endl;
            std::cout << " Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
            std::cout << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
            ++nb_holes;
        }
    }
    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;

    std::ofstream out("filled.off");
    out.precision(17);
    out << poly << std::endl;

    JMeshOFFImporter mimp;
    JMeshPtr newmesh = mimp.readFile("filled.off");
    return newmesh;
}


