#ifndef FUNFEM_UTILHPP
#define FUNFEM_UTILHPP

#include "expression.hpp"


template <typename Mesh>
std::vector<size_t> getBoundaryDof( const CutFESpace<Mesh> &Vh, const TimeSlab &In, std::list<int> label = {}) {

    bool all_label = (label.size() == 0);
    std::vector<size_t> dof2set;
    dof2set.clear();
    const auto &cutTh(Vh.get_mesh());

    // Loop through all fitted boundary elements 
    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
        int k                 = idxK[0];

        assert(idxK.size() == 1);
        
        const auto &K(cutTh.Th[kb]);
        const auto &BE(cutTh.be(idx_be));
        const auto &FK(Vh[k]);

        if (util::contain(label, BE.lab) || all_label) {
            // for( int ic=0; ic<Vh.N;++ic) {
            
            // DOF for BDM is only one component (v*n)
            for (int ic = 0; ic < 1; ++ic) {

                // Loop over degrees of freedom of current element
                for (int df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df) {

                    int id_item = FK.DFOnWhat(df);

                    if (id_item < K.nv) {
                        assert(0);
                    } else if (id_item < K.nv + K.ne) {
                        // std::cout << " on edge  " <<FK.DFOnWhat(df) << std::endl;
                        int id_face = id_item - K.nv;
                        if (id_face == ifac) {
                            for (int dof_time=0; dof_time < In.NbDoF(); ++dof_time) {
                                size_t df_glob = FK.loc2glb(df, dof_time);
                                dof2set.push_back(df_glob);
                            }
                            
                        }
                    } else {
                        // std::cout << " on face  " << FK.DFOnWhat(df) << std::endl;
                    }
                }
            }
        }
    }

    return dof2set;
}


template <typename Mesh>
void setBoundaryDof( const std::vector<size_t>& dof2set, const FunFEM<Mesh> &gh,std::span<double> b) {
for (const auto& df : dof2set)
    b[df] = gh(df);
}

inline void setBoundaryDof( const std::vector<size_t>& dof2set, double val ,std::span<double> b) {
for (const auto& df : dof2set)
    b[df] = val;
}

inline void setBoundaryDof(const size_t N, const std::vector<size_t>& dof2set, double val ,std::map<std::pair<int,int>,double>& A) {

    // auto [N,M] = size(A);
    std::map<std::pair<int,int>,double> C;
    std::map<std::pair<int, int>, double> P;
    for (int i = 0; i < N; ++i) {
        P[std::make_pair(i, i)] = 1;
    }

    for (auto &i0 : dof2set) {
        P[std::make_pair(i0, i0)] = 0;
    }

    SparseMatrixRC<double> AA(N, N, A);
    SparseMatrixRC<double> PP(N, N, P);
    multiply(PP, AA, C);
    A = std::move(C);
    // SparseMatrixRC<double> CC(N,N,C);
    // multiply(CC, PP, A);
    

    for (auto &i0 : dof2set) {
        A[std::make_pair(i0, i0)] = val;
    }

}


// /**
//  * Strong boundary conditions for BDM1 in time
// */
// template <typename Mesh>
// void setBoundaryDof(const FunFEM<Mesh> &gh, const ActiveMesh<Mesh> &cutTh, const TimeSlab &In, std::span<double> b, std::list<int> label) {

//     assert(b.size() >= gh.size());

//     bool all_label = (label.size() == 0);
//     // std::map<int, double> dof2set;
//     const auto &Vh(gh.getSpace());

//     // Loop through all fitted boundary elements 
//     for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
//          idx_be += cutTh.next_boundary_element()) {

//         int ifac;
//         const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
//         std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
//         int k                 = idxK[0];

//         assert(idxK.size() == 1);
        
//         const auto &K(cutTh.Th[kb]);
//         const auto &BE(cutTh.be(idx_be));
//         const auto &FK(Vh[k]);

//         if (util::contain(label, BE.lab) || all_label) {

//             // for( int ic=0; ic<Vh.N;++ic) {
            
//             // DOF for BDM is only one component (v*n)
//             for (int ic = 0; ic < 1; ++ic) {

//                 // Loop over degrees of freedom of current element
//                 for (int df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df) {

//                     int id_item = FK.DFOnWhat(df);

//                     if (id_item < K.nv) {
//                         assert(0);
//                     } else if (id_item < K.nv + K.ne) {
//                         // std::cout << " on edge  " <<FK.DFOnWhat(df) << std::endl;
//                         int id_face = id_item - K.nv;
//                         if (id_face == ifac) {

//                             for (int dof_time; dof_time < In.NbDoF(); ++dof_time) {
//                                 int df_glob = FK.loc2glb(df, dof_time);
//                                 // dof2set.insert({df_glob, gh(df_glob)});
//                                 // dof2set.insert({df_glob, gh(df_glob)});

//                                 b[df_glob] = gh(df_glob);

//                                 // std::cout << df_glob << std::endl;    
//                             }
                            
//                         }
//                     } else {
//                         // std::cout << " on face  " << FK.DFOnWhat(df) << std::endl;
//                     }
//                 }
//             }
//         }
//         // getchar();
//     }

// //     assert(this->pmat_.size() == 1);
// //     eraseAndSetRow(this->get_nb_dof(), *(this->pmat_[0]), this->rhs_, dof2set);
// }

#endif