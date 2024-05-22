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
#ifndef PARALLEL_DUMMY_CFMPI_HPP
#define PARALLEL_DUMMY_CFMPI_HPP

#include "../num/util.hpp"

// Dummy mpi class for seq compilation
class MPIcf {
  public:
    MPIcf(int &argc, char **&argv) {}
    MPIcf() {}

    static int my_rank() { return 0; }
    static int size() { return 1; }

    static bool IamMaster() { return true; }
    static int Master() { return 0; }

    static inline double Wtime() { return CPUtime(); };

    static inline const int first_element(int n) { return 0; }
    static inline const int next_element(int n) { return 1; }
    static inline const int last_element(int n) { return n; }

    template <typename T> static inline void Bcast(T &a, int who, int size) {}
};
#endif // PARALLEL_DUMMY_CFMPI_HPP
