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

#ifndef FORCE_PROJECT_UTIL_HPP
#define FORCE_PROJECT_UTIL_HPP

namespace force_project {

template <typename T>
bool checkIfValid(double xc, double yc, double delta_l, const std::vector<std::shared_ptr<T>> &levelSets_v) {

    for (int i = 0; i < levelSets_v.size(); i++) {
        double val = levelSets_v[i]->operator()(R2(xc, yc));
        if (val < delta_l)
            return false;
    }
    return true;
}

} // namespace force_project

#endif