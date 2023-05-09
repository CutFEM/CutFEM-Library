
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

#ifndef COMMON_BASE_TIME_INTERFACE_HPP_
#define COMMON_BASE_TIME_INTERFACE_HPP_

#include <iostream>
#include <cassert>
#include <bitset>
#include <memory>
#include "../concept/function.hpp"
#include "interface_levelSet.hpp"
#include "AlgoimInterface.hpp"

template <typename Mesh> class TimeInterface {
  public:
    using interface_t = Interface<Mesh>;
    using mesh_t      = Mesh;

  private:
    std::vector<std::unique_ptr<interface_t>> interface_;
    int n_;
    const QuadratureFormular1d *time_quadrature_;

  public:
    TimeInterface(const QuadratureFormular1d &qTime) : interface_(qTime.n), n_(qTime.n), time_quadrature_(&qTime) {}

    TimeInterface(const QuadratureFormular1d *qTime) : interface_(qTime->n), n_(qTime->n), time_quadrature_(qTime) {}

    TimeInterface(int nt) : interface_(nt), n_(nt), time_quadrature_(Lobatto(exactLobatto_nPt(nt))) {}

    /// @brief Copy constructor is removed
    TimeInterface(const TimeInterface &) = delete;

    /// @brief Assigment is removed
    void operator=(const TimeInterface &) = delete;

    /// @brief Move constructor
    TimeInterface(TimeInterface &&v) = default;

    /// @brief Move assignment
    TimeInterface &operator=(TimeInterface &&v) = default;

    /// @brief Destructor
    ~TimeInterface() = default;

    template <typename Fct> void init(int i, const Mesh &Th, const Fct &ls) {
        assert(0 <= i && i < n_);
        if constexpr (typeFunFEM<Fct>) {
            interface_[i] = std::make_unique<InterfaceLevelSet<mesh_t>>(Th, ls);
        } else {
            interface_[i] = std::make_unique<AlgoimInterface<MeshQuad2, Fct>>(Th, ls);
        }
    }
    // template <typename Fct> void init(const Mesh &Th, const KN<Fct> &ls);
    // {
    //     assert(n_ == ls.size());
    //     for (int i = 0; i < n_; ++i) {
    //         //         interface_[i] = std::make_unique<interface_t>(Th, ls[i]);
    //         if constexpr (typeFunFEM<Fct>) {
    //             interface_[i] = std::make_unique < InterfaceLevelset<mesh_t>(Th, ls[i]);
    //         } else {
    //             interface_[i] = std::make_unique < AlgoimInterface<MeshQuad2, Fct>(Th, ls[i]);
    //         }
    //     }
    // }

    interface_t *operator[](int i) const {
        assert(0 <= i && i < n_);
        return interface_[i].get();
    }
    interface_t *operator()(int i) const {
        assert(0 <= i && i < n_);
        return interface_[i].get();
    }

    int size() const { return n_; }
    const QuadratureFormular1d *get_quadrature_time() const { return time_quadrature_; }

    const std::vector<std::unique_ptr<interface_t>> &interface() const { return interface_; }
};
#endif