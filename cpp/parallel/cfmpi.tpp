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

template <typename T> inline long MPIcf::Bcast(T &a, int who, int size) { return WBcast(&a, size, who, Communicator_); }

template <typename T> inline long MPIcf::Bcast(T &a, int who, int size, const MPI_Comm &comm) {
    return WBcast(&a, size, who, comm);
}

template <typename T> inline long MPIcf::Bcast(const KN<T> &a, int who) {
    assert(&a);
    CheckContigueKN(a);
    return WBcast((T *)a, a.N(), who, Communicator_);
}
template <typename T> inline long MPIcf::Bcast(std::span<T> a, int who) {
    assert(a.data());
    return WBcast(a.data(), a.size(), who, Communicator_);
}
template <typename T> inline long MPIcf::Bcast(const KN<T> &a, int who, const MPI_Comm &comm) {
    assert(&a);
    CheckContigueKN(a);
    return WBcast((T *)a, a.N(), who, comm);
}

template <typename T> inline void MPIcf::Reduce(const T &myData, T &globalData, int size, const MPI_Op &op, int root) {
    MPI_Reduce(&myData, &globalData, size, MPI_TYPE<T>::TYPE(), op, root, Communicator_);
}

template <typename T> inline void MPIcf::Reduce(const KN<T> &myData, KN<T> &globalData, const MPI_Op &op, int root) {
    MPI_Reduce((T *)(myData), (T *)globalData, myData.size(), MPI_TYPE<T>::TYPE(), op, root, Communicator_);
}

template <typename T>
inline void MPIcf::Reduce(std::span<T> myData, std::span<T> globalData, const MPI_Op &op, int root) {
    MPI_Reduce(myData.data(), globalData.data(), myData.size(), MPI_TYPE<T>::TYPE(), op, root, Communicator_);
}

template <typename T> inline void MPIcf::AllReduce(const T &myData, T &globalData, const MPI_Op &op) {
    MPI_Allreduce(&myData, &globalData, 1, MPI_TYPE<T>::TYPE(), op, Communicator_);
}

template <typename T> inline void MPIcf::AllReduce(const T &myData, T &globalData, int size, const MPI_Op &op) {
    MPI_Allreduce(&myData, &globalData, size, MPI_TYPE<T>::TYPE(), op, Communicator_);
}

template <typename T> inline void MPIcf::AllReduce(const T *myData, T *globalData, int size, const MPI_Op &op) {
    MPI_Allreduce(myData, globalData, size, MPI_TYPE<T>::TYPE(), op, Communicator_);
}

template <typename T> inline void MPIcf::AllReduce(const KN<T> &myData, KN<T> &globalData, const MPI_Op &op) {
    CheckContigueKN(myData);
    MPI_Allreduce((T *)(myData), globalData, myData.size(), MPI_TYPE<T>::TYPE(), op, Communicator_);
}

template <typename T> inline void MPIcf::Scan(const T &myData, T &globalData, const MPI_Op &op) {
    MPI_Scan(&myData, &globalData, 1, MPI_TYPE<T>::TYPE(), op, Communicator_);
}
