/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef COMMON_SORT_ARRAY_HPP
#define COMMON_SORT_ARRAY_HPP

template <typename T, int N> struct SortArray {};

template <typename T> struct SortArray<T, 1> {
   T v[1];
   SortArray(T *a) { v[0] = a[0]; }
   SortArray(const T &a0) { v[0] = a0; }
   SortArray() {}
   bool operator==(const SortArray<T, 1> &t) const { return v[0] == t.v[0]; }
   bool operator<(const SortArray<T, 1> &t) const { return v[0] < t.v[0]; }
   size_t hash() const { return (size_t)v[0]; }
};

template <typename T> struct SortArray<T, 2> {
   //  using std::std::swap;
   T v[2];
   SortArray(T *a) {
      v[0] = a[0];
      v[1] = a[1];
      if (v[0] > v[1])
         std::swap(v[0], v[1]);
   }
   SortArray(const T &a0, const T &a1) {
      v[0] = a0;
      v[1] = a1;
      if (v[0] > v[1])
         std::swap(v[0], v[1]);
   }
   SortArray() {}
   bool operator==(const SortArray<T, 2> &t) const {
      return v[0] == t.v[0] && v[1] == t.v[1];
   }
   bool operator<(const SortArray<T, 2> &t) const {
      return v[0] != t.v[0] ? v[0] < t.v[0] : v[1] < t.v[1];
   }
   T operator[](const int i) const {
      assert(i < 2);
      return v[i];
   }
   size_t hash() const { return (size_t)v[0]; }
};

template <typename T> struct SortArray<T, 3> {
   T v[3];
   SortArray(T *a) {
      v[0] = a[0];
      v[1] = a[1];
      v[2] = a[2];
      if (v[0] > v[1])
         std::swap(v[0], v[1]);
      if (v[1] > v[2]) {
         std::swap(v[1], v[2]);
         if (v[0] > v[1])
            std::swap(v[0], v[1]);
         assert(v[0] <= v[1] && v[1] <= v[2]);
      }
   }

   SortArray(const T &a0, const T &a1, const T &a2) {
      v[0] = a0;
      v[1] = a1;
      v[2] = a2;
      if (v[0] > v[1])
         std::swap(v[0], v[1]);
      if (v[1] > v[2]) {
         std::swap(v[1], v[2]);
         if (v[0] > v[1])
            std::swap(v[0], v[1]);
         assert(v[0] <= v[1] && v[1] <= v[2]);
      }
   }

   SortArray() {}
   bool operator==(const SortArray<T, 3> &t) const {
      return v[0] == t.v[0] && v[1] == t.v[1] && v[2] == t.v[2];
   }

   bool operator<(const SortArray<T, 3> &t) const {
      return v[0] != t.v[0] ? v[0] < t.v[0]
                            : (v[1] != t.v[1] ? v[1] < t.v[1] : v[2] < t.v[2]);
   }
   T operator[](const int i) const {
      assert(i < 3);
      return v[i];
   }
   size_t hash() const { return (size_t)v[0]; }
};

template <typename T> struct SortArray<T, 4> {
   T v[4];
   SortArray(T *a) {
      v[0] = a[0];
      v[1] = a[1];
      v[2] = a[2];
      v[3] = a[3];
      if (v[0] > v[1])
         std::swap(v[0], v[1]);
      if (v[1] > v[2]) {
         std::swap(v[1], v[2]);
         if (v[0] > v[1])
            std::swap(v[0], v[1]);
      }
      if (v[2] > v[3]) {
         std::swap(v[2], v[3]);
         if (v[1] > v[2]) {
            std::swap(v[1], v[2]);
            if (v[0] > v[1])
               std::swap(v[0], v[1]);
         }
         assert(v[0] <= v[1] && v[1] <= v[2] && v[2] <= v[3]);
      }
   }

   SortArray(const T &a0, const T &a1, const T &a2, const T &a3) {
      v[0] = a0;
      v[1] = a1;
      v[2] = a2;
      v[3] = a3;
      if (v[0] > v[1])
         std::swap(v[0], v[1]);
      if (v[1] > v[2]) {
         std::swap(v[1], v[2]);
         if (v[0] > v[1])
            std::swap(v[0], v[1]);
      }
      if (v[2] > v[3]) {
         std::swap(v[2], v[3]);
         if (v[1] > v[2]) {
            std::swap(v[1], v[2]);
            if (v[0] > v[1])
               std::swap(v[0], v[1]);
         }
         assert(v[0] <= v[1] && v[1] <= v[2] && v[2] <= v[3]);
      }
   }

   SortArray() {}
   bool operator==(const SortArray<T, 4> &t) const {
      return v[0] == t.v[0] && v[1] == t.v[1] && v[2] == t.v[2] &&
             v[3] == t.v[3];
   }

   bool operator<(const SortArray<T, 4> &t) const {
      return v[0] != t.v[0]
                 ? v[0] < t.v[0]
                 : (v[1] != t.v[1]
                        ? v[1] < t.v[1]
                        : ((v[2] != t.v[2]) ? v[2] < t.v[2] : v[3] < t.v[3]));
   }
   T operator[](const int i) const {
      assert(i < 4);
      return v[i];
   }
   size_t hash() const { return (size_t)v[0]; }
};

template <typename T, int N>
bool operator==(const SortArray<T, N> &v1, const SortArray<T, N> &v2) {
   for (int j = 0; j < N; j++) {
      if (v1.v[j] != v2.v[j]) {
         return false;
      }
   }
   return true;
}

template <typename T, int N>
std::ostream &operator<<(std::ostream &f, const SortArray<T, N> &item) {
   for (int i = 0; i < N; ++i)
      f << " " << item.v[i];
   return f;
}

#endif
