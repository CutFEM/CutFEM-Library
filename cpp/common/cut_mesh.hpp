/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/

#ifndef COMMON_CUT_MESH_HPP
#define COMMON_CUT_MESH_HPP

#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "base_interface.hpp"
#include "time_interface.hpp"

template <typename E> struct Cut_Part {
    static const int dim = E::RdHat::d;               ///< dimension of the element
    typedef SortArray<Ubyte, dim + 1> ElementIdx;     ///< the vertices of a triangle of the cut:
    typedef const ElementIdx *const_element_iterator; ///< const iterator for the elements
    typedef ElementIdx *element_iterator;             ///< iterator for the elements
    typedef typename E::Rd Rd;                        ///< type for the physical coordinates
    typedef typename E::RdHat RdHat;                  ///< type for the coordinates in the reference element

    const Virtual_Partition<E> *partition_; ///< the partition
    const Partition<E> ip;                  ///< the partition in the reference element
    const Physical_Partition<E> pp;         ///< the partition in the physical element
    int sign_cut_;                          ///< the sign of the cut

    Cut_Part(const Partition<E> p, int s);          ///< constructor with a reference partition
    Cut_Part(const Physical_Partition<E> p, int s); ///< constructor with a physical partition
    Cut_Part(const Cut_Part<E> &p);                 ///< copy constructor

    Cut_Part &operator=(const Cut_Part<E> &p); // copy operator

    // GETTERS
    int get_sign() const;                                        ///< return the sign of the cut
    int get_sign_node(int i) const;                              ///< return the sign of the i-th node
    int get_nb_element() const;                                  ///< return the number of elements in the cut
    int get_local_domain_id() const;                             ///< return the local domain id
    void get_list_node(std::vector<typename E::Rd> &node) const; ///< return the list of nodes
    Rd get_vertex(const_element_iterator it, const int i) const; ///< return the i-th vertex of an element
    CutElement<E> get_element(int k) const;                      ///< return the k-th element

    // OTHER METHODS
    // GIVE THE MEASURE OF THE CUT PART IN Rd
    double measure() const;                          ///< return the measure of the cut part
    double measure(const_element_iterator it) const; ///< return the measure of an element

    // //GIVE THE MEASURE OF CUT PART OF A FACE IN RdBord
    // double measureBord(int ifac) const {return
    // partition_->measureBord(sign_cut_, ifac);}

    bool multi_interface() const; ///< return true if there are more than one interface

    Rd mapToPhysicalElement(const_element_iterator it,
                            const RdHat Phat) const; ///< map a point in the reference element to the physical element

    // ITERATORS
    const_element_iterator element_begin() const; ///< return an iterator to the beginning of the elements
    const_element_iterator element_end() const;   ///< return an iterator to the end of the elements
    const_element_iterator other_side_element_begin()
        const; ///< return an iterator to the beginning of the elements on the other side of the cut
    const_element_iterator
    other_side_element_end() const; ///< return an iterator to the end of the elements on the other side of the cut
};

/**
 * @brief
 *
 * @note The active mesh of a given subdomain corresponds to those elements of the
 * background mesh that have a non-empty intersection with the subdomain. In the
 * space-time method, the active mesh is taken as all elements that intersect the
 * subdomain during any time instance in the time-slab In=[t_{n-1}, t_{n}]. This
 * class allows for several types of active meshes, such as multiple active meshes
 * of different subdomains, single active mesh of a single subdomain, or active meshes
 * of surfaces.
 * @tparam Mesh Mesh
 */
template <typename Mesh> class ActiveMesh {

    // Public Types
  public:
    typedef typename Mesh::Element Element;             ///< Mesh element
    typedef typename Element::Face Face;                ///< Mesh face
    typedef typename Mesh::Rd Rd;                       ///< Physical coordinate tuple
    typedef typename Mesh::BorderElement BorderElement; ///< Mesh border element
    typedef typename Element::RdHat RdHat;              ///< For parametrization
    typedef SortArray<int, 2> pairIndex;                //?
    static const int nea = Element::nea; ///<  Number of adjacent elements (4 in Tet,  3 in Tria, 2 in seg)
    static const int nva = Element::nva; ///<  Number of vertex in adjacent face (3 in Tet, 2 in Tria, 1 in seg)

    const Mesh &Th; ///< background mesh

    // Public Attributes
  public:
    /**
     * @brief Given index of subdomain and index of element in active mesh, get corresponding index in background mesh.
     * @param subdomain subdomain index
     * @param idxK_ActiveMesh index of element in active mesh
     * @param idxK_backMesh index of element in background mesh
     * @return Calling \code{.cpp} idx_in_background_mesh_[subdomain](idxK_ActiveMesh)\endcode returns
     * \code{.cpp}idxK_backMesh\endcode
     */
    std::vector<std::vector<int>> idx_in_background_mesh_;

    /**
     * @brief Given index of subdomain and index of element in background mesh, get index in
     * element in the active mesh corresponding to the subdomain
     * @param subdomain subdomain index
     * @param idxK_ActiveMesh index of element in active mesh
     * @param idxK_backMesh index of element in background mesh
     * @return Calling \code{.cpp} idx_from_background_mesh_[subdomain](idxK_backMesh)\endcode returns
     * \code{.cpp}idxK_ActiveMesh\endcode
     */
    std::vector<std::map<int, int>> idx_from_background_mesh_;

    /**
     * @brief Vector of the interfaces in the time quadrature points itq //?
     * @note interface_id_[itq] returns what exactly ?
     */
    std::vector<std::map<std::pair<int, int>, std::vector<std::pair<const Interface<Mesh> *, int>>>>
        interface_id_;                   // [time_quad](domain_id, idx_k) ->
                                         // [n_interface](interface, sign)
    std::vector<int> idx_element_domain; //! What does this do?

    // number of quadrature points to approximate integral over the time slab In
    int nb_quadrature_time_{1};

    // a list of the elements that are in the active mesh for some time quadrature point, but not for all
    std::vector<std::vector<std::map<int, bool>>> in_active_mesh_; // [dom][itq][idx_element] -> true/false

    // Constructors
  public:
    /**
     * @brief Construct an ActiveMesh object of the whole background mesh.
     * @note Useful constructor for adding subdomains.
     * @param th Background mesh.
     */
    ActiveMesh(const Mesh &th);

    /**
     * @brief Construct two active meshes, one for each subdomain on the positive and negative side of the level set.
     * @param th Background mesh.
     * @param interface Interface object.
     */
    ActiveMesh(const Mesh &th, const Interface<Mesh> &interface);

    ActiveMesh(const Mesh &th, const TimeInterface<Mesh> &interface);

    void truncate(const Interface<Mesh> &interface, int sign_domain);

    /**
     * @brief Create active mesh for one time-dependent subdomain.
     * @note When creating one bulk domain, one removes (truncates) the other domain
     * corresponding to the other sign of the signed distance function.
     * @param interface Time-dependent interface.
     * @param sign_domain Sign of the signed distance function in the domain you want to remove.
     */
    void truncate(const TimeInterface<Mesh> &interface, int sign_domain);

    /**
     * @brief I don't know what this does. //?
     *
     * @param interface
     * @param sign_domain
     */
    void add(const Interface<Mesh> &interface, int sign_domain);

    /**
     * @brief Create an active mesh of a stationary surface.
     *
     * @param interface Interface<Mesh> object.
     */
    void createSurfaceMesh(const Interface<Mesh> &interface);

    /**
     * @brief Create an active mesh of a time-dependent surface.
     *
     * @param interface TimeInterface<Mesh> object.
     */
    void createSurfaceMesh(const TimeInterface<Mesh> &interface);

    /**
     * @brief Dont know what this does. //?
     *
     * @param interface
     */
    void addArtificialInterface(const Interface<Mesh> &interface);

  public:
    int nbElmts() const;
    int NbElement() const;
    int nbBrdElmts() const;
    int nbVertices() const;

    /**
     * @brief Get the Element of a given index in the active mesh
     *
     * @tparam Mesh
     * @param i index in the active mesh
     * @return const ActiveMesh<Mesh>::Element&
     */
    const Element &operator[](int i) const;
    const BorderElement &be(int i) const;

    int get_nb_domain() const;
    int get_nb_element(int i) const;
    int get_nb_element() const;
    int get_domain_element(const int k) const;

    bool isCut(int k, int t) const;
    bool isCutFace(int k, int ifac, int t) const;
    bool isStabilizeElement(int k) const;
    bool isInactive(int k, int t) const;

    const Interface<Mesh> &get_interface(int k, int t) const;

    Partition<Element> get_partition(int k, int t) const;

    int get_sign_cut(int k, int t) const;

    Cut_Part<Element> get_cut_part(int k, int t) const;

    Cut_Part<typename Element::Face> get_cut_face(Face &face, int k, int ifac, int t) const;

    std::vector<int> idxAllElementFromBackMesh(int k, int d) const;

    std::vector<int> getAllDomainId(int k) const;

    int idxElementFromBackMesh(int k) const;
    int idxElementFromBackMesh(int k, int i) const;
    int idxElementInBackMesh(const int k) const;
    int idxElementInBackMesh(const int k, int i) const;
    int ElementAdj(const int k, int &j) const;
    std::tuple<int, int> elementAdjacent(const int k, const int i) const;

    void info() const;

    virtual int first_element() const { return MPIcf::first_element(this->get_nb_element()); }
    virtual int next_element() const { return MPIcf::next_element(this->get_nb_element()); }
    virtual int last_element() const { return MPIcf::last_element(this->get_nb_element()); }

    virtual int first_boundary_element() const { return MPIcf::my_rank(); }
    virtual int next_boundary_element() const { return MPIcf::size(); }
    virtual int last_boundary_element() const { return this->Th.nbBrdElmts(); }


    // Private Member Functions
  private:
    void init(const Interface<Mesh> &interface);
    void init(const TimeInterface<Mesh> &interface);
    bool check_exist(int k, int dom) const;
    int idxK_begin(int i) const;
    int idxK_in_domain(int k, int i) const;
    Physical_Partition<Element> build_local_partition(const int k, int t = 0) const;
    Physical_Partition<Face> build_local_partition(Face &face, const int k, int ifac, int t = 0) const;
};

typedef ActiveMesh<Mesh2> ActiveMeshT2;
typedef ActiveMesh<Mesh3> ActiveMeshT3;
typedef ActiveMesh<MeshQuad2> ActiveMeshQ2;
typedef ActiveMesh<MeshHexa> ActiveMeshQ3;

#include "cut_mesh.tpp"

#endif
