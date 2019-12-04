#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>


//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;

//namespace SMS = CGAL::Surface_mesh_simplification ;

//int main( int argc, char** argv )
//{
//  if (argc<3)
//  {
//    std::cerr << "Usage: " << argv[0] << " input.off minimal_edge_length [out.off]\n";
//    return EXIT_FAILURE;
//  }

//  Surface_mesh surface_mesh;

//  std::ifstream is(argv[1]) ; is >> surface_mesh ;
//  double threshold = atof(argv[2]);

//  if (!CGAL::is_triangle_mesh(surface_mesh)){
//    std::cerr << "Input geometry is not triangulated." << std::endl;
//    return EXIT_FAILURE;
//  }

//  int r = SMS::edge_collapse
//            (surface_mesh
//             , CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double>(threshold)
//             , CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh))
//                               .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface_mesh))
//                               .get_cost (SMS::Edge_length_cost <Surface_mesh>())
//                               .get_placement(SMS::Midpoint_placement<Surface_mesh>())
//            );

//  std::cout << "\nFinished...\n" << r << " edges removed.\n"
//            << (surface_mesh.size_of_halfedges()/2) << " final edges.\n" ;

//  std::ofstream os( argc > 3 ? argv[3] : "out.off" ) ; os << surface_mesh ;

//  return EXIT_SUCCESS ;
//}
template <class Refs, class T, class Norm>
class MEPP_Common_Facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
    protected:
        // tag
        int m_tag;

        // normal
        Norm m_normal;

        // color
        float m_color[3];

    public:
        // life cycle
        MEPP_Common_Facet()
        {
            color(0.5f, 0.5f, 0.5f);
        }

        // tag
        const int& tag() const { return m_tag;  }
        int& tag() { return m_tag;  }
        void tag(const int& t)  { m_tag = t; }

        // normal
        typedef Norm Normal_3;
        Normal_3& normal() { return m_normal; }
        const Normal_3& normal() const { return m_normal; }

        // color
        float color(int index) { return m_color[index]; }
        void color(float r, float g, float b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MEPP_Common_Halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
    protected:
        // tag
        int m_tag;

#if (0)
        // normal
        Norm m_normal;  // AJOUT C�line : halfedge normal (sum of 2 incident facet normals)
#endif

        // option for edge superimposing
        bool m_control_edge;

        // texture coordinates : AJOUT Laurent Chevalier
        float m_texture_coordinates[2];

    public:
        // life cycle
        MEPP_Common_Halfedge()
        {
            m_control_edge = true;

            // texture coordinates : AJOUT Laurent Chevalier
            texture_coordinates(0.0f, 0.0f);
        }

        // tag
        const int& tag() const { return m_tag;  }
        int& tag() { return m_tag;  }
        void tag(const int& t)  { m_tag = t; }

#if (0)
        // normal : AJOUT C�line
        typedef Norm Normal_3;
        Normal_3& normal() { return m_normal; }
        const Normal_3& normal() const { return m_normal; }
#endif

        // texture coordinates : AJOUT Laurent Chevalier
        float texture_coordinates(int index) { return m_texture_coordinates[index]; }
        void texture_coordinates(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }

        // control edge
        bool& control_edge()  { return m_control_edge; }
        const bool& control_edge()  const { return m_control_edge; }
        void control_edge(const bool& flag) { m_control_edge  = flag; }
};

template <class Refs, class T, class P, class Norm>
class MEPP_Common_Vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
    protected:
        // tag
        int m_tag;

        // normal
        Norm m_normal;

        // color
        float m_color[3];

        // texture coordinates
        float m_texture_coordinates[2];

    public:
        // life cycle
        MEPP_Common_Vertex()
        {
            color(0.5f, 0.5f, 0.5f);
            texture_coordinates(0.0f, 0.0f);
        }
        // repeat mandatory constructors
        MEPP_Common_Vertex(const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
        {
            color(0.5f, 0.5f, 0.5f);
            texture_coordinates(0.0f, 0.0f);
        }

        // color
        float color(int index) { return m_color[index]; }
        void color(float r, float g, float b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }

        // texture coordinates
        float texture_coordinates(int index) { return m_texture_coordinates[index]; }
        void texture_coordinates(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }

        // normal
        typedef Norm Normal_3;
        //typedef Norm Vector;

        Normal_3& normal() { return m_normal; }
        const Normal_3& normal() const { return m_normal; }

        // tag
        int& tag() {  return m_tag; }
        const int& tag() const {  return m_tag; }
        void tag(const int& t)  { m_tag = t; }
};
template <class Refs, class T, class P, class Norm, class Plane>
class Enriched_facet :
    /*************** HERITAGE FACETTE ***************/
//	#include <polyhedron_enrichment_facet.h>
    /*************************************************/
    virtual public MEPP_Common_Facet<Refs, T, Norm>
{
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Enriched_halfedge :
    /*************** HERITAGE HALFEDGE ***************/
//	#include <polyhedron_enrichment_halfedge.h>
    /*************************************************/
    virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
};

// a refined vertex with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Enriched_vertex :
    /*************** HERITAGE VERTEX ***************/
//	#include <polyhedron_enrichment_vertex.h>
    /*************************************************/
    virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
    public:
        Enriched_vertex()
        {
        }

        // Le constructeur PAS par d�faut appelle les constructeurs PAR DEFAUT des classes ancetres
        // en plus de celles appelles appell�es explicitement.
        // On a besoin d'appeller explicitement le constructeur(pt) de base pour la cr�ation de polyhedre
        Enriched_vertex(const P& pt) : MEPP_Common_Vertex<Refs, T, P, Norm>(pt)
        {
            this->point() = pt;
        }

        // La creation du polyhedre implique un appel au constructeur par copie,
        // qui appelle le constructeur par d�faut de base et ne copie pas le point.
        // Il faut donc avoir un constructeur par copie explicite qui s'occupe du point.
        Enriched_vertex(const Enriched_vertex& v)
        {
            this->point() = v.point();
        }
};

struct Enriched_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template <class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_vertex<Refs,
                          CGAL::Tag_true,
                          Point,
                          Normal> Vertex;
    };

    // wrap face
    template <class Refs, class Traits>
    struct Face_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef typename Traits::Plane_3 Plane;
        typedef Enriched_facet<Refs,
                         CGAL::Tag_true,
                         Point,
                         Normal,
                         Plane> Face;
    };

    // wrap halfedge
    template <class Refs, class Traits>
    struct Halfedge_wrapper
    {
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_halfedge<Refs,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            Normal> Halfedge;
    };
};

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel,Enriched_items> Surface;
int main(int argc, char *argv[])
{
    if (argc<3)
      {
        std::cerr << "Usage: " << argv[0] << " input.off minimal_edge_length [out.off]\n";
        return EXIT_FAILURE;
      }
    Surface surface_mesh;
    CGAL::Surface_mesh_simplification::Count_stop_predicate<Surface> stop(std::stoi(argv[2]));

    std::ifstream is(argv[1]) ; is >> surface_mesh ;
    Surface * surface = &surface_mesh;
    int r = CGAL::Surface_mesh_simplification::edge_collapse
                (*surface
                ,stop,
                 CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, *surface))
                                 .halfedge_index_map(get(CGAL::halfedge_external_index, *surface)));
    std::cout << "\nFinished...\n" << r << " edges removed.\n"
                << (surface_mesh.size_of_halfedges()/2) << " final edges.\n" ;

      std::ofstream os( argc >4  ? argv[3] : "out.off" ) ; os << surface_mesh ;


    return 0;
}
