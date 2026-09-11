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
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef CUTFEM_AGGREGATE_EXTENSION_HPP
#define CUTFEM_AGGREGATE_EXTENSION_HPP

#include "expression.hpp"
#include "macroElement.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace cutfem {

/** Mark fully physical, uncut cells as roots for data known only in the physical domain. */
template <typename Mesh>
std::vector<bool> fullyInteriorRootEligibility(const ActiveMesh<Mesh> &Th, int time_index = 0) {
    std::vector<bool> eligible(Th.get_nb_element(), false);
    for (int k = 0; k < Th.get_nb_element(); ++k) {
        eligible[k] = !Th.isInactive(k, time_index) && !Th.isCut(k, time_index);
    }
    return eligible;
}

/** Physical volume fractions obtained from the active mesh's geometric cut partition. */
template <typename Mesh>
std::vector<double> cutPartPhysicalVolumeFractions(const ActiveMesh<Mesh> &Th, int time_index = 0) {
    std::vector<double> fraction(Th.get_nb_element(), 0.);
    for (int k = 0; k < Th.get_nb_element(); ++k) {
        if (Th.isInactive(k, time_index)) {
            continue;
        }
        if (!Th.isCut(k, time_index)) {
            fraction[k] = 1.;
            continue;
        }

        const Cut_Part<typename Mesh::Element> physical_part(Th.get_cut_part(k, time_index));
        fraction[k] = std::clamp(physical_part.measure() / Th[k].measure(), 0., 1.);
    }
    return fraction;
}

/**
 * Stationary aggregate partition for discrete extensions from a physical
 * domain to an active mesh.
 *
 * A cell is a root when it is eligible and its caller-provided physical volume
 * fraction is at least root_fraction. Every other active cell is assigned to
 * the nearest root in the same domain by a face-neighbour breadth-first
 * search. Supplying root eligibility separately from the volume fraction lets
 * a geometry implementation exclude cut cells even if numerical quadrature
 * rounds their volume fraction to one.
 */
template <typename Mesh> class StationaryAggregateMacro : public GMacro {
  public:
    StationaryAggregateMacro(const ActiveMesh<Mesh> &Th, const std::vector<double> &physical_fraction,
                             double root_fraction, std::vector<bool> root_eligible = {})
        : Th_(Th), physical_fraction_(physical_fraction), root_fraction_(root_fraction),
          root_eligible_(std::move(root_eligible)) {
        if (physical_fraction_.size() != static_cast<std::size_t>(Th_.get_nb_element())) {
            throw std::invalid_argument(
                "StationaryAggregateMacro: one volume fraction per active element is required");
        }
        if (!(0. < root_fraction_ && root_fraction_ <= 1.)) {
            throw std::invalid_argument("StationaryAggregateMacro: root fraction must lie in (0,1]");
        }
        if (root_eligible_.empty()) {
            root_eligible_.assign(physical_fraction_.size(), true);
        } else if (root_eligible_.size() != physical_fraction_.size()) {
            throw std::invalid_argument(
                "StationaryAggregateMacro: one root-eligibility flag per active element is required");
        }
        createRootMap();
    }

    double physicalFraction(int k) const { return physical_fraction_.at(k); }
    std::size_t numberOfRootCells() const { return number_of_root_cells_; }

    int maximumChainLength() const {
        int maximum = 0;
        for (const auto &entry : small_element) {
            maximum = std::max(maximum, entry.second.chain_position);
        }
        return maximum;
    }

  private:
    const ActiveMesh<Mesh> &Th_;
    std::vector<double> physical_fraction_;
    double root_fraction_;
    std::vector<bool> root_eligible_;
    std::size_t number_of_root_cells_ = 0;

    void createRootMap() {
        const int number_of_elements = Th_.get_nb_element();
        std::vector<int> root(number_of_elements, -1);
        std::vector<int> distance(number_of_elements, -1);
        std::queue<int> frontier;

        for (int k = 0; k < number_of_elements; ++k) {
            const double fraction = physical_fraction_[k];
            if (!std::isfinite(fraction) || fraction < -1e-12 || fraction > 1. + 1e-10) {
                throw std::runtime_error("StationaryAggregateMacro: invalid physical volume fraction");
            }

            if (root_eligible_[k] && fraction >= root_fraction_) {
                root[k] = k;
                distance[k] = 0;
                frontier.push(k);
                ++number_of_root_cells_;
            } else {
                small_element[k] = SmallElement(k);
                small_element[k].area = std::max(0., fraction) * Th_[k].measure();
            }
        }

        if (frontier.empty()) {
            throw std::runtime_error("StationaryAggregateMacro: the active mesh contains no eligible root cell");
        }

        while (!frontier.empty()) {
            const int k = frontier.front();
            frontier.pop();

            for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) {
                int neighbour_face = ifac;
                const int kn = Th_.ElementAdj(k, neighbour_face);
                if (kn < 0 || distance[kn] >= 0) {
                    continue;
                }
                if (Th_.get_domain_element(kn) != Th_.get_domain_element(k)) {
                    continue;
                }

                root[kn] = root[k];
                distance[kn] = distance[k] + 1;
                frontier.push(kn);
            }
        }

        for (auto &[k, cell] : small_element) {
            if (root[k] < 0) {
                throw std::runtime_error(
                    "StationaryAggregateMacro: an active component is disconnected from every eligible root cell");
            }
            cell.setRoot(root[k]);
            cell.setChainPosition(distance[k]);

            auto insertion = macro_element.try_emplace(
                root[k], root[k], physical_fraction_[root[k]] * Th_[root[k]].measure(), this);
            insertion.first->second.add(k, cell.area);
        }

        // std::map insertions do not invalidate references. Populate this
        // lookup only after all macro elements have been created.
        for (auto &entry : macro_element) {
            auto &aggregate = entry.second;
            for (const int k : aggregate.idx_element) {
                idx_element_to_macro_element[k] = &aggregate;
            }
        }
    }
};

struct AggregateExtensionInfo {
    std::size_t root_dofs = 0;
    std::size_t extended_dofs = 0;
};

struct AggregateNodalError {
    double absolute = 0.;
    double reference_scale = 0.;

    double relative() const { return absolute / std::max(1., reference_scale); }
};

namespace detail {

template <typename Mesh>
void validateAggregateExtensionSpace(std::span<double> coefficients, const CutFESpace<Mesh> &Wh) {
    if (coefficients.size() != static_cast<std::size_t>(Wh.get_nb_dof())) {
        throw std::invalid_argument("aggregate extension: coefficient vector has the wrong size");
    }

    std::size_t cell_local_dofs = 0;
    for (int k = 0; k < Wh.get_nb_element(); ++k) {
        cell_local_dofs += Wh[k].NbDoF();
    }
    if (cell_local_dofs != coefficients.size()) {
        throw std::invalid_argument("aggregate extension: the extension space must be discontinuous");
    }
}

} // namespace detail

/**
 * Fill small-cell coefficients by canonical polynomial extension from each
 * aggregate root:
 *
 *     w_K(x_a) = w_R(x_a).
 *
 * Root coefficients must already be initialized. Wh must be a nodal,
 * discontinuous space with one independently writable set of coefficients per
 * active cell. The polynomial degree and number of components are otherwise
 * selected by the caller.
 */
template <typename Mesh>
std::size_t extendAggregateCoefficients(std::span<double> coefficients, const CutFESpace<Mesh> &Wh,
                                        const GMacro &aggregates) {
    detail::validateAggregateExtensionSpace(coefficients, Wh);

    std::size_t extended_dofs = 0;
    for (const auto &[k, cell] : aggregates.small_element) {
        const auto &FK = Wh[k];
        const auto &FR = Wh[cell.index_root];

        if (FK.N != FR.N || FK.NbDoF() != FR.NbDoF()) {
            throw std::runtime_error("extendAggregateCoefficients: nonuniform extension space");
        }

        KNMK<double> root_basis(FR.NbDoF(), FR.N, 1);
        for (int component = 0; component < FK.N; ++component) {
            const int target_begin = FK.dfcbegin(component);
            const int target_end = FK.dfcend(component);
            const int root_begin = FR.dfcbegin(component);
            const int root_end = FR.dfcend(component);
            const int interpolation_points = FK.tfe->NbPtforInterpolation;

            if (target_end - target_begin != interpolation_points || root_end - root_begin != interpolation_points) {
                throw std::runtime_error(
                    "extendAggregateCoefficients: expected one nodal degree of freedom per interpolation point");
            }

            for (int i = target_begin; i < target_end; ++i) {
                const int interpolation_point = i - target_begin;
                const auto physical_point = FK.Pt(interpolation_point);
                const auto root_reference_point = FR.T.mapToReferenceElement(physical_point);
                FR.BF(Fop_D0, root_reference_point, root_basis);

                double value = 0.;
                for (int j = root_begin; j < root_end; ++j) {
                    value += coefficients[FR.loc2glb(j)] * root_basis(j, component, op_id);
                }
                coefficients[FK.loc2glb(i)] = value;
                ++extended_dofs;
            }
        }
    }

    return extended_dofs;
}

/** Interpolate pointwise data on aggregate roots without sampling non-roots. */
template <typename Mesh, typename Function>
std::size_t interpolateAggregateRootCoefficients(std::span<double> coefficients, const CutFESpace<Mesh> &Wh,
                                                 const GMacro &aggregates, Function &&source) {
    detail::validateAggregateExtensionSpace(coefficients, Wh);

    std::size_t root_dofs = 0;
    auto &&source_ref = source;
    for (int k = 0; k < Wh.get_nb_element(); ++k) {
        if (aggregates.isSmall(k)) {
            continue;
        }

        const auto &FK = Wh[k];
        const int domain = Wh.get_domain(k);
        const int number_of_interpolation_points = FK.tfe->NbPtforInterpolation;
        KNM<double> values(number_of_interpolation_points, FK.N);
        KN<double> local_coefficients(FK.NbDoF());

        for (int interpolation_point = 0; interpolation_point < number_of_interpolation_points;
             ++interpolation_point) {
            auto physical_point = FK.Pt(interpolation_point);
            for (int component = 0; component < FK.N; ++component) {
                values(interpolation_point, component) =
                    std::invoke(source_ref, physical_point, component, domain);
            }
        }

        FK.Pi_h(values, local_coefficients);
        for (int i = 0; i < FK.NbDoF(); ++i) {
            coefficients[FK.loc2glb(i)] = local_coefficients(i);
            ++root_dofs;
        }
    }
    return root_dofs;
}

/** Interpolate on aggregate roots and extend the resulting polynomials. */
template <typename Mesh, typename Function>
AggregateExtensionInfo interpolateAggregateExtension(std::span<double> coefficients, const CutFESpace<Mesh> &Wh,
                                                      const GMacro &aggregates, Function &&source) {
    AggregateExtensionInfo info;
    info.root_dofs = interpolateAggregateRootCoefficients(
        coefficients, Wh, aggregates, std::forward<Function>(source));
    info.extended_dofs = extendAggregateCoefficients(coefficients, Wh, aggregates);
    return info;
}

/**
 * Own the aggregate map, coefficients, and discrete function produced by one
 * stationary physical-to-active extension. The finite element space remains
 * caller-selected; this keeps the operation independent of a particular PDE.
 */
template <typename Mesh> class AggregateExtension {
  public:
    template <typename Function>
    AggregateExtension(const CutFESpace<Mesh> &Wh, const std::vector<double> &physical_fraction,
                       double root_fraction, std::vector<bool> root_eligible, Function &&source)
        : Wh_(Wh), aggregates_(Wh.get_mesh(), physical_fraction, root_fraction, std::move(root_eligible)),
          coefficients_(Wh.get_nb_dof(), 0.), function_(Wh, coefficients_) {
        info_ = interpolateAggregateExtension(std::span<double>(coefficients_), Wh_, aggregates_,
                                              std::forward<Function>(source));
    }

    AggregateExtension(const AggregateExtension &) = delete;
    AggregateExtension &operator=(const AggregateExtension &) = delete;
    AggregateExtension(AggregateExtension &&) = delete;
    AggregateExtension &operator=(AggregateExtension &&) = delete;

    FunFEM<Mesh> &function() { return function_; }
    const FunFEM<Mesh> &function() const { return function_; }
    const StationaryAggregateMacro<Mesh> &aggregates() const { return aggregates_; }
    const AggregateExtensionInfo &info() const { return info_; }

    double constraintResidual() const {
        double residual = 0.;
        for (const auto &[k, small_cell] : aggregates_.small_element) {
            const auto &FK = Wh_[k];
            for (int interpolation_point = 0; interpolation_point < FK.tfe->NbPtforInterpolation;
                 ++interpolation_point) {
                const auto physical_point = FK.Pt(interpolation_point);
                for (int component = 0; component < FK.N; ++component) {
                    residual = std::max(
                        residual,
                        std::fabs(function_.eval(k, physical_point, component, op_id) -
                                  function_.eval(small_cell.index_root, physical_point, component, op_id)));
                }
            }
        }
        return residual;
    }

    void writeDiagnostics(std::ostream &output, std::string_view label = "Aggregate extension") const {
        output << label << ": roots=" << aggregates_.numberOfRootCells()
               << ", aggregates=" << aggregates_.macro_element.size()
               << ", extended_cells=" << aggregates_.small_element.size()
               << ", extended_dofs=" << info_.extended_dofs
               << ", max_chain=" << aggregates_.maximumChainLength()
               << ", constraint_residual=" << constraintResidual() << '\n';
    }

    template <typename Function> AggregateNodalError nodalError(Function &&reference) const {
        AggregateNodalError error;
        auto &&reference_function = reference;
        for (int k = 0; k < Wh_.get_nb_element(); ++k) {
            const auto &FK = Wh_[k];
            const int domain = Wh_.get_domain(k);
            for (int interpolation_point = 0; interpolation_point < FK.tfe->NbPtforInterpolation;
                 ++interpolation_point) {
                auto physical_point = FK.Pt(interpolation_point);
                for (int component = 0; component < FK.N; ++component) {
                    const double exact_value =
                        std::invoke(reference_function, physical_point, component, domain);
                    error.reference_scale = std::max(error.reference_scale, std::fabs(exact_value));
                    error.absolute = std::max(
                        error.absolute,
                        std::fabs(function_.eval(k, physical_point, component, op_id) - exact_value));
                }
            }
        }
        return error;
    }

  private:
    const CutFESpace<Mesh> &Wh_;
    StationaryAggregateMacro<Mesh> aggregates_;
    std::vector<double> coefficients_;
    FunFEM<Mesh> function_;
    AggregateExtensionInfo info_;
};

/** Construct an aggregate extension using the active mesh's geometric Cut_Part. */
template <typename Mesh, typename Function>
AggregateExtension<Mesh> makeCutPartAggregateExtension(const CutFESpace<Mesh> &Wh, Function &&source,
                                                       double root_fraction = 1., int time_index = 0) {
    return AggregateExtension<Mesh>(
        Wh, cutPartPhysicalVolumeFractions<Mesh>(Wh.get_mesh(), time_index), root_fraction,
        fullyInteriorRootEligibility<Mesh>(Wh.get_mesh(), time_index), std::forward<Function>(source));
}

} // namespace cutfem

#endif
