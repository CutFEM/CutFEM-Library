/**
 * @file cpp/problem/time_scheme.hpp
 * @brief
 *
 * @copyright

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

#ifndef CUTFEM_TIME_SCHEME_HPP
#define CUTFEM_TIME_SCHEME_HPP
namespace Scheme {

namespace RK3 {

template <class InputIt, class OutputIt>
void step1(InputIt first0, InputIt last0, InputIt first1, double dt, OutputIt r_first) {
    while (first0 != last0) {
        *r_first++ = *first0++ + dt * *first1++;
    }
}

template <class InputIt, class OutputIt>
void step2(InputIt first0, InputIt last0, InputIt first1, InputIt first2, double dt, OutputIt r_first) {
    while (first0 != last0) {
        *r_first++ = 3. / 4 * *first0++ + 1. / 4 * *first1++ + 1. / 4 * dt * *first2++;
    }
}

template <class InputIt, class OutputIt>
void step3(InputIt first0, InputIt last0, InputIt first1, InputIt first2, double dt, OutputIt r_first) {
    while (first0 != last0) {
        *r_first++ = 1. / 3 * *first0++ + 2. / 3 * *first1++ + 2. / 3 * dt * *first2++;
    }
}
} // namespace RK3
} // namespace Scheme

#endif