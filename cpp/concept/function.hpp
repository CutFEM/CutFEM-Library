/**
 * @file cpp/concept/function.hpp
 *
 * @brief Contains concepts related to functions
 *
 */

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

#ifndef CONCEPT_FUNCTION_HPP
#define CONCEPT_FUNCTION_HPP

#include <type_traits>
#include <concepts>
#include <span>
#include <vector>
#include <functional>

class Mesh1;
class Mesh2;
class Mesh3;
class MeshQuad2;
class MeshHexa;
template <typename C>
concept typeMesh = std::is_same_v<Mesh1, C> || std::is_same_v<Mesh2, C> || std::is_same_v<Mesh3, C> ||
                   std::is_same_v<MeshHexa, C> || std::is_same_v<MeshQuad2, C>;

template <typename C>
concept meshQuadrilateral = std::is_same_v<MeshQuad2, C> || std::is_same_v<MeshHexa, C>;

// template <typename L>
// concept AlgoimLevelset = requires(L phi) {
//     phi<typename T>(const algoim::uvector<T, typename N> &);

// }

template <typename M> class GFESpace;
template <typename M> class CutFESpace;
template <typename C>
concept Space = std::is_same_v<GFESpace<typename C::Mesh>, C> || std::is_same_v<CutFESpace<typename C::Mesh>, C>;

// Concepts for functions
// template <int D, int E, typename Vec, typename T> class Function;
template <typename M> class FunFEM;
template <typename C>
concept typeFunFEM = std::is_same_v<FunFEM<typename C::Mesh>, C>;

class R1;
class R2;
class R3;
template <typename C>
concept Node = std::is_same_v<R1, C> || std::is_same_v<R2, C> || std::is_same_v<R3, C>;

using fct_scalar_ptr = std::add_pointer_t<double(double *)>;
using fct_scalar_R2  = std::add_pointer_t<double(R2)>;
using fct_scalar_R3  = std::add_pointer_t<double(R3)>;
template <typename fct_t>
concept FunctionScalar =
    std::is_same_v<fct_t, fct_scalar_ptr> || std::is_same_v<fct_t, fct_scalar_R2> ||
    std::is_same_v<fct_t, fct_scalar_R3> || std::is_convertible_v<fct_t, std::function<double(std::span<double>)>>;

using fct_ptr_int = std::add_pointer_t<double(double *, int)>;
using fct_R2_int  = std::add_pointer_t<double(R2, int)>;
using fct_R3_int  = std::add_pointer_t<double(R3, int)>;

template <typename fct_t>
concept FunctionLevelSet =
    std::is_same_v<fct_t, fct_ptr_int> || std::is_same_v<fct_t, fct_R2_int> || std::is_same_v<fct_t, fct_R3_int> ||
    std::is_convertible_v<fct_t, std::function<double(std::span<double>, int)>>;

using fct_ptr_int_double = std::add_pointer_t<double(double *, int, const double)>;
using fct_R2_int_double  = std::add_pointer_t<double(R2, int, const double)>;
using fct_R3_int_double  = std::add_pointer_t<double(R3, int, double)>;

template <typename fct_t>
concept FunctionLevelSetTime = std::is_same_v<fct_t, fct_ptr_int_double> || std::is_same_v<fct_t, fct_R2_int_double> ||
                               std::is_same_v<fct_t, fct_R3_int_double>;

using fct_ptr_int_int = std::add_pointer_t<double(double *, int, int)>;
using fct_R2_int_int  = std::add_pointer_t<double(R2, int, int)>;
using fct_R3_int_int  = std::add_pointer_t<double(R3, int, int)>;

template <typename fct_t>
concept FunctionDomain = std::is_same_v<fct_t, fct_ptr_int_int> || std::is_same_v<fct_t, fct_R2_int_int> ||
                         std::is_same_v<fct_t, fct_R3_int_int> ||
                         std::is_convertible_v<fct_t, std::function<double(std::span<double>, int, int)>>;

using fct_ptr_int_int_double = std::add_pointer_t<double(double *, int, int, double)>;
using fct_R2_int_int_double  = std::add_pointer_t<double(R2, int, int, double)>;
using fct_R3_int_int_double  = std::add_pointer_t<double(R3, int, int, double)>;

template <typename fct_t>
concept FunctionDomainTime =
    std::is_same_v<fct_t, fct_ptr_int_int_double> || std::is_same_v<fct_t, fct_R2_int_int_double> ||
    std::is_same_v<fct_t, fct_R3_int_int_double>;


// using fct_ptr_int_double = std::add_pointer_t<double(double *, int, double)>;
// using fct_R2_int_double  = std::add_pointer_t<double(R2, int, double)>;
// using fct_R3_int_double  = std::add_pointer_t<double(R3, int, double)>;

// template <typename fct_t>
// concept FunctionTime =
//     std::is_same_v<fct_t, fct_ptr_int_double> || std::is_same_v<fct_t, fct_R2_int_double> ||
//     std::is_same_v<fct_t, fct_R3_int_double>;

template <typename R> class KN_;
template <typename R> class KN;
template <typename C> struct is_vector : std::false_type {};
template <typename T> struct is_vector<std::vector<T>> : std::true_type {};
template <typename T> struct is_vector<std::span<T>> : std::true_type {};
template <typename T> struct is_vector<KN<T>> : std::true_type {};
template <typename T> struct is_vector<KN_<T>> : std::true_type {};
template <typename C>
concept Vector = is_vector<C>::value;

template <typename C>
concept NonAllocVector =
    std::is_same_v<C, KN_<typename C::element_type>> || std::is_same_v<C, std::span<typename C::element_type>>;

// // Concept for shared_ptr
// template <typename T> struct is_shared_ptr : std::false_type {};
// template <typename T>
// struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};
// template <typename T>
// concept IsSharedPtr = is_shared_ptr<T>::value;

// // Concepts for expression
// template <typename T> class Expression;
// template <typename C>
// concept IsExpression = std::derived_from<C, Expression<typename C::T>>;

#endif