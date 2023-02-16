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

class Mesh1;
class Mesh2;
class Mesh3;
class MeshQuad2;
class MeshHexa;
template <typename C>
concept typeMesh = std::is_same_v<Mesh1, C> || std::is_same_v<Mesh2, C> ||
                   std::is_same_v<Mesh3, C> || std::is_same_v<MeshHexa, C> ||
                   std::is_same_v<MeshQuad2, C>;

// Concepts for functions
// template <int D, int E, typename Vec, typename T> class Function;
template <typename M> class FunFEM;
template <typename C>
concept typeFunFEM = std::is_same_v<FunFEM<typename C::Mesh>, C>;

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