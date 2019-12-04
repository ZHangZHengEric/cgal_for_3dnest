#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

namespace VSA = CGAL::Surface_mesh_approximation;
namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
    // read input surface triangle mesh
    Mesh mesh;
    std::ifstream file("data/join_ok.off");
    file >> mesh;
    std::cout<<"12123"<<std::endl;
    // output indexed triangles
    std::vector<Kernel::Point_3> anchors;
    std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices

    // free function interface with named parameters
    bool is_manifold = VSA::approximate_triangle_mesh(mesh,
                                                      CGAL::parameters::seeding_method(VSA::HIERARCHICAL). // hierarchical seeding
                                                      verbose_level(CGAL::Surface_mesh_approximation::VERBOSE).
                                                      max_number_of_proxies(3000). // seeding with maximum number of proxies
                                                      number_of_iterations(50). // set number of clustering iterations after seeding

                                                      anchors(std::back_inserter(anchors)). // anchor vertices
                                                      triangles(std::back_inserter(triangles))); // indexed triangles

    std::cout << "#anchor vertices: " << anchors.size() << std::endl;
    std::cout << "#triangles: " << triangles.size() << std::endl;

    std::ofstream output("output.off");

    output<<"OFF"<<std::endl;
    output<<anchors.size()<<" "<<triangles.size()<<" "<<0<<std::endl;
    for(int i=0; i< anchors.size() ;i++){
        output<<anchors[i].x()<<" "<<anchors[i].y()<<" "<<anchors[i].z()<<std::endl;
    }
    for(int i =0; i< triangles.size() ;i++){
        output<< 3<<" "<<triangles[i][0]<<" "<<triangles[i][1]<<" "<<triangles[i][2]<<std::endl;
    }
    output.close();
    if (is_manifold) {
        std::cout << "oriented, 2-manifold output." << std::endl;




        // convert from soup to surface mesh
        PMP::orient_polygon_soup(anchors, triangles);
        Mesh output;
        PMP::polygon_soup_to_polygon_mesh(anchors, triangles, output);
        if (CGAL::is_closed(output) && (!PMP::is_outward_oriented(output)))
            PMP::reverse_face_orientations(output);

        std::ofstream out("dump.off");
        out << output;
        out.close();
    }

    return EXIT_SUCCESS;
}
