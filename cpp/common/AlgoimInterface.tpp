
#ifndef ALGOIM_INTERFACE_TPP
#define ALGOIM_INTERFACE_TPP

template <typeMesh M, typename L>
AlgoimInterface<M, L>::AlgoimInterface(const M &Mesh, L &phi_, int label) : Interface<M>(Mesh), phi(phi_) {
    std::cout << "AlgoimInterface constructor"
              << "\n";
    make_algoim_patch(label);
}

template <typeMesh M, typename L> void AlgoimInterface<M, L>::make_algoim_patch(int label) {

    std::cout << "make_algoim_patch  method"
              << "\n";
    using mesh_t  = M;
    using Element = typename AlgoimInterface<M, L>::Element;
    assert(this->backMesh);
    this->faces_.resize(0); // reinitialize arrays
    this->vertices_.resize(0);
    this->element_of_face_.resize(0);
    this->outward_normal_.resize(0);
    this->face_of_element_.clear();

    const mesh_t &Th = *(this->backMesh); // background mesh

    // Iterate over all elements in the background mesh
    for (int k = 0; k < Th.nbElmts(); k++) {

        const Element &K(Th[k]);

        // Get coordinates of current quadrilateral
        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

        if (q.nodes.size() == 0) {
            // K is not cut
            continue;
        }

        else {
            // K is cut
            cut_elements.insert({k, q});
            number_of_cut_elements += 1;
        }
    }
}

template <typeMesh M, typename L>
SignElement<typename AlgoimInterface<M, L>::Element> AlgoimInterface<M, L>::get_SignElement(int k) const {
    using Element = typename AlgoimInterface<M, L>::Element;

    double loc_ls[Element::nv];
    for (int i = 0; i < Element::nv; ++i) {
        loc_ls[i] = phi(this->backMesh->operator[](k).at(i));
    }
    return SignElement<Element>(loc_ls);
}

template <typeMesh M, typename L>
Partition<typename AlgoimInterface<M, L>::Element> AlgoimInterface<M, L>::get_partition(int k) const {
    using Element = typename AlgoimInterface<M, L>::Element;

    //assert(0);
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

template <typeMesh M, typename L> double AlgoimInterface<M, L>::measure(const AlgoimInterface<M,L>::Face &f) const {
    assert(0);
    return 0.;
}

template <typeMesh M, typename L>
typename AlgoimInterface<M, L>::Rd
AlgoimInterface<M,L>::mapToPhysicalFace(int ifac, const typename AlgoimInterface<M,L>::Element::RdHatBord x) const {
    //typename AlgoimInterface<M,L>::Rd N[nve];
    assert(0);
//     for (int i = 0; i < nve; ++i)
//         N[i] = this->vertices_[this->faces_[ifac][i]];
//     return geometry::map_point_to_simplex(N, x);
    return typename AlgoimInterface<M, L>::Rd();
}

// if index is in the cut_elements map, its corresponding element is cut
template <typeMesh M, typename L>
bool AlgoimInterface<M,L>::isCut(int k) const {
    return (cut_elements.find(k) != cut_elements.end());
}
#endif