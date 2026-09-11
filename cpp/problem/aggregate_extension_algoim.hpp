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

#ifndef CUTFEM_AGGREGATE_EXTENSION_ALGOIM_HPP
#define CUTFEM_AGGREGATE_EXTENSION_ALGOIM_HPP

#include "../FESpace/aggregate_extension.hpp"
#include "baseProblem.hpp"

#ifdef USE_LAPACK

namespace cutfem {

/** Physical volume fractions obtained with the caller's Algoim volume rule. */
template <typename Mesh, typename LevelSet, typename Option>
std::vector<double> algoimPhysicalVolumeFractions(const ActiveMesh<Mesh> &Th, LevelSet &phi, const Option &option,
                                                  int time_index = 0, double time = 0.) {
    std::vector<double> fraction(Th.get_nb_element(), 0.);
    phi.setTime(time);
    for (int k = 0; k < Th.get_nb_element(); ++k) {
        if (Th.isInactive(k, time_index)) {
            continue;
        }
        if (!Th.isCut(k, time_index)) {
            fraction[k] = 1.;
            continue;
        }

        const int domain = Th.get_domain_element(k);
        phi.setElementFromBackMesh(Th.idxElementInBackMesh(k), domain);
        const auto quadrature = quadGenVol(Th[k], phi, option, domain);
        double physical_measure = 0.;
        for (const double weight : quadrature.weights) {
            physical_measure += weight;
        }
        fraction[k] = std::clamp(physical_measure / Th[k].measure(), 0., 1.);
    }
    return fraction;
}

/** Construct an aggregate extension using high-order Algoim physical volume fractions. */
template <typename Mesh, typename LevelSet, typename Option, typename Function>
AggregateExtension<Mesh> makeAlgoimAggregateExtension(const CutFESpace<Mesh> &Wh, LevelSet &phi,
                                                       const Option &option, Function &&source,
                                                       double root_fraction = 1., int time_index = 0,
                                                       double time = 0.) {
    return AggregateExtension<Mesh>(
        Wh, algoimPhysicalVolumeFractions<Mesh>(Wh.get_mesh(), phi, option, time_index, time), root_fraction,
        fullyInteriorRootEligibility<Mesh>(Wh.get_mesh(), time_index), std::forward<Function>(source));
}

} // namespace cutfem

#endif

#endif
