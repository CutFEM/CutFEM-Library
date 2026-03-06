
#ifndef ALGOIM_INTERFACE_TPP
#define ALGOIM_INTERFACE_TPP

template <typeMesh M, typename L>
AlgoimInterface<M, L>::AlgoimInterface(const M &Mesh, const L &phi_, int label) : Interface<M>(Mesh), phi(phi_) {

    make_algoim_patch(label);
}

template <typeMesh M, typename L> void AlgoimInterface<M, L>::make_algoim_patch(int label) {

    using mesh_t  = M;
    using Element = typename AlgoimInterface<M, L>::Element;
    assert(this->backMesh);
    this->faces_.resize(1); // reinitialize arrays
    this->vertices_.resize(0);
    this->element_of_face_.resize(0);
    // this->outward_normal_.resize(0);
    this->face_of_element_.clear();

    ProblemOption options;
    options.order_space_element_quadrature_ = 3;
    options.algoim_bernstein_deg_ = 2;

    const mesh_t &Th = *(this->backMesh); // background mesh

    // Iterate over all elements in the background mesh
    for (int k = 0; k < Th.nbElmts(); k++) {

        const Element &K(Th[k]);

        // // Get coordinates of current quadrilateral
        // const auto &V0(K.at(0)); // vertex 0
        // const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

        // algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        // algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        // algoim::QuadratureRule<2> q =
        //     algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

        if constexpr (std::is_same_v<std::remove_cvref_t<L>, FunFEM<M>>) {
            phi.setElementFromBackMesh(k, 0);
        }
        auto q = quadGenSurf(K, phi, options);
        // double cut_threshold = 1e-6;


        //if ((q.points.size() == 0) || (std::accumulate(q.weights.begin(), q.weights.end(), 0.0) < cut_threshold)) {
        if (q.size() == 0) {
            // K is not cut
            continue;
        } else {
            // K is cut
            this->element_of_face_.push_back(k);    // add element to list of cut elements
            cut_elements.insert({k, q});
            number_of_cut_elements += 1;
        }
    }
}

/**
 * @brief Finds the sign of the interface in the nodes of the element with index k
 * 
 * @tparam M Mesh type
 * @tparam L Algoim level set type
 * @param k index of element in background mesh
 * @return SignElement<typename AlgoimInterface<M, L>::Element> Array of sign of level set in the nodes
 */
template <typeMesh M, typename L>
SignElement<typename AlgoimInterface<M, L>::Element> AlgoimInterface<M, L>::get_SignElement(int k) const {
    using Element = typename AlgoimInterface<M, L>::Element;

    if constexpr (std::is_same_v<std::remove_cvref_t<L>, FunFEM<M>>) {
        phi.setElement(k);
    }
    double loc_ls[Element::nv]; // list with size = number of vertices on the element
    for (int i = 0; i < Element::nv; ++i) {
        loc_ls[i] = phi(this->backMesh->operator[](k).at(i));
    }
    return SignElement<Element>(loc_ls);
}

template <typeMesh M, typename L>
Partition<typename AlgoimInterface<M, L>::Element> AlgoimInterface<M, L>::get_partition(int k) const {
    using Element = typename AlgoimInterface<M, L>::Element;

    if constexpr (std::is_same_v<std::remove_cvref_t<L>, FunFEM<M>>) {
        phi.setElement(k);
    }
    // assert(0);
    double loc_ls[Element::nv];
    for (int i = 0; i < Element::nv; ++i) {
        loc_ls[i] = phi(this->backMesh->operator[](k).at(i));
    }

    return Partition<Element>((*this->backMesh)[k], loc_ls);
}

template <typeMesh M, typename L>
Partition<typename AlgoimInterface<M, L>::Element::Face>
AlgoimInterface<M, L>::get_partition_face(const typename AlgoimInterface<M, L>::Element::Face &face, int k,
                                          int ifac) const {
    using Element = typename AlgoimInterface<M, L>::Element;

    assert(0);
    double loc_ls[Element::Face::nv];
    // for (int i = 0; i < Element::Face::nv; ++i) {
    //     int j     = Element::nvhyperFace[ifac][i];
    //     int iglb  = this->backMesh->at(k, j);
    //     loc_ls[i] = ls_[iglb];
    // }
    return Partition<typename Element::Face>(face, loc_ls);
}

template <typeMesh M, typename L>
void AlgoimInterface<M, L>::cut_partition(Physical_Partition<typename AlgoimInterface<M, L>::Element> &local_partition,
                                          std::vector<ElementIdx> &new_element_idx, std::list<int> &erased_element,
                                          int sign_part) const {
    std::cout << " An element might be cut multiplue time, and it is not "
                 "suppose to happen"
              << std::endl;
    exit(EXIT_FAILURE);
};

template <typeMesh M, typename L> double AlgoimInterface<M, L>::measure(const AlgoimInterface<M, L>::Face &f) const {
    assert(0);
    return 0.;
}

template <typeMesh M, typename L> double AlgoimInterface<M, L>::measure(int i) const {
    assert(0);
    return 0.;
}

template <typeMesh M, typename L>
typename AlgoimInterface<M, L>::Rd AlgoimInterface<M, L>::normal(int k, std::span<double> x) const {
    assert(0);
    if constexpr (std::is_same_v<std::remove_cvref_t<L>, FunFEM<M>>) {
        phi.setElement(k);
    }
    return phi.normal(x);
}

template <typeMesh M, typename L>
typename AlgoimInterface<M, L>::Rd
AlgoimInterface<M, L>::mapToPhysicalFace(int ifac, const typename AlgoimInterface<M, L>::Element::RdHatBord x) const {
    // typename AlgoimInterface<M,L>::Rd N[nve];
    assert(0);
    //     for (int i = 0; i < nve; ++i)
    //         N[i] = this->vertices_[this->faces_[ifac][i]];
    //     return geometry::map_point_to_simplex(N, x);
    return typename AlgoimInterface<M, L>::Rd();
}

// if index is in the cut_elements map, its corresponding element is cut
template <typeMesh M, typename L> bool AlgoimInterface<M, L>::isCut(int k) const {
    return (cut_elements.find(k) != cut_elements.end());
}


//! REWRITE FOR GENERAL MESHES
template <typeMesh M, typename L> bool AlgoimInterface<M, L>::isCutFace(int k, int ifac) const {

    const M &Th = *(this->backMesh); // background mesh

    const Element &K(Th[k]);

    // If the face is on the (background) mesh boundary, consider it not an interior cut-face
    // (no shared neighbor). ElementAdj returns -1 when adjacent element is outside the mesh.
    // Note: ElementAdj modifies ifac by reference, so we use a copy
    int ifac_copy = ifac;
    int kn = Th.ElementAdj(k, ifac_copy);
    if (kn == -1) {
        std::cout << "Face " << ifac << " of bg element " << k << " is on the boundary of the mesh, not a cut face." << std::endl;
        return false;
    }

    ProblemOption options;
    options.algoim_surface_quad_deg_ = 3 ;
    options.algoim_bernstein_deg_    = 2;
    
    if constexpr (std::is_same_v<std::remove_cvref_t<L>, FunFEM<M>>) {
        phi.setElement(k);
    }
    const auto quad_rule = quadGenFace(K, phi, options, ifac);

    // // Get coordinates of current quadrilateral
    // const auto &V0(K.at(0)); // vertex 0
    // const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

    // algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    // algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    // // Map face index to (dim, side) for quadGen.
    // // Based on GenericElement<DataQuad2>::nvedge = {{0,1},{1,2},{2,3},{3,0}}
    // // and vertex coordinates R2::QuadHat = { (0,0),(1,0),(1,1),(0,1) }
    // // the edges correspond to:
    // //  ifac==0 : edge {0,1} -> bottom  (y == ymin) -> dim=1, side=0
    // //  ifac==1 : edge {1,2} -> right   (x == xmax) -> dim=0, side=1
    // //  ifac==2 : edge {2,3} -> top     (y == ymax) -> dim=1, side=1
    // //  ifac==3 : edge {3,0} -> left    (x == xmin) -> dim=0, side=0
    // int dim = -1, side = -1;
    // switch (ifac) {
    //     case 0: dim = 1; side = 0; break; // bottom
    //     case 1: dim = 0; side = 1; break; // right
    //     case 2: dim = 1; side = 1; break; // top
    //     case 3: dim = 0; side = 0; break; // left
    //     default: assert(false && "Unexpected face index in isCutFace");
    // }

    // assert(0);
    // return false;
    // // algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), dim, side,
    // //                                                   quadrature_order);

    // // return (q.nodes.size() != 0);
    const bool is_cut_face = (quad_rule.size() != 0);
    
    // std::cout << "Face " << ifac << " of bg element " << k << (is_cut_face ? " is cut." : " is not cut.") << std::endl;

    return is_cut_face;
}


#endif