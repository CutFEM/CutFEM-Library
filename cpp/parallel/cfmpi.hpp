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
#ifndef _CF_MPI_HPP_
#define _CF_MPI_HPP_

#include "/opt/homebrew/Cellar/open-mpi/4.1.4_2/include/mpi.h"
#include "../num/util.hpp"
#include "../common/RNM.hpp"
#include <cassert>

// get the MPI datatype
template <class T> struct MPI_TYPE {
   static MPI_Datatype TYPE() { return MPI_BYTE; }
};
;
template <> struct MPI_TYPE<long> {
   static MPI_Datatype TYPE() { return MPI_LONG; }
};
template <> struct MPI_TYPE<int> {
   static MPI_Datatype TYPE() { return MPI_INT; }
};
template <> struct MPI_TYPE<Uint> {
   static MPI_Datatype TYPE() { return MPI_UNSIGNED; }
};
template <> struct MPI_TYPE<Ulint> {
   static MPI_Datatype TYPE() { return MPI_UNSIGNED_LONG; }
};
template <> struct MPI_TYPE<Usint> {
   static MPI_Datatype TYPE() { return MPI_UNSIGNED_SHORT; }
};
template <> struct MPI_TYPE<double> {
   static MPI_Datatype TYPE() { return MPI_DOUBLE; }
};
template <> struct MPI_TYPE<char> {
   static MPI_Datatype TYPE() { return MPI_CHAR; }
};
template <> struct MPI_TYPE<byte> {
   static MPI_Datatype TYPE() { return MPI_CHAR; }
};

// The Tag for the messages
template <class T> struct MPI_TAG {};
template <> struct MPI_TAG<long> {
   static const int TAG = 5;
};
template <> struct MPI_TAG<double> {
   static const int TAG = 4;
};
template <> struct MPI_TAG<KN<long> *> {
   static const int TAG = 11;
};
template <> struct MPI_TAG<KN<double> *> {
   static const int TAG = 12;
};
template <> struct MPI_TAG<KNM<long> *> {
   static const int TAG = 14;
};
template <> struct MPI_TAG<KNM<double> *> {
   static const int TAG = 15;
};

// for syncro communication
static MPI_Request *Syncro_block = reinterpret_cast<MPI_Request *>(1);
const size_t sizempibuf          = 1024 * 320;

template <class T> void CheckContigueKNM(const KNM_<T> &t) {
   if (t.step != 1 && !t.IsVector1()) {
      std::cout << " step= " << t.step << " size " << t.N() << " " << &t[0]
                << " " << &t[1] << std::endl;
      assert(0 && "Sorry the array is not contigue (step != 1) ");
   }
}

template <class T> void CheckContigueKN(const KN_<T> &t) {
   if (t.step != 1 && t.N() > 1) {
      std::cout << " step= " << t.step << " size " << t.N() << " " << &t[0]
                << " " << &t[1] << std::endl;
      assert(0 && "Sorry the array is not contigue (step != 1) ");
   }
}

template <class R>
long WSend(R *v, int l, int who, int tag, MPI_Comm comm, MPI_Request *rq) {
   long ret = 0;
   MPI_Request rq0, *request = &rq0;

   if (rq == Syncro_block || rq == 0)
      ret = MPI_Send((void *)v, l, MPI_TYPE<R>::TYPE(), who, tag, comm);
   else {
      ret =
          MPI_Isend((void *)v, l, MPI_TYPE<R>::TYPE(), who, tag, comm, request);
      if (rq)
         *rq = *request;
      else
         MPI_Request_free(request);
   }
   return ret;
}

template <class R>
long WRecv(R *v, int n, int who, int tag, MPI_Comm comm, MPI_Request *rq) {
   MPI_Status status;
   if (rq && (rq != Syncro_block))
      return MPI_Irecv(reinterpret_cast<void *>(v), n, MPI_TYPE<R>::TYPE(), who,
                       tag, comm, rq);
   else
      return MPI_Recv(reinterpret_cast<void *>(v), n, MPI_TYPE<R>::TYPE(), who,
                      tag, comm, &status);
}

template <class R> long WBcast(R *v, int n, int who, MPI_Comm comm) {
   assert(v && n > 0);
   return MPI_Bcast(reinterpret_cast<void *>(v), n, MPI_TYPE<R>::TYPE(), who,
                    comm);
}

// manage group of processors
class MPIcf {
 private:
   static int my_rank_;
   static int size_;
   static const MPI_Comm &Communicator_;
   static MuteStdOstream *mute_;
   static MPI_Request *rq;
   static bool usePetsc_;
   static int loopDivision;

 public:
   MPIcf(int &argc, char **&argv);
   MPIcf();

   static int my_rank() { return my_rank_; }
   static int size() { return size_; }

   static bool IamMaster() { return my_rank_ == 0; }
   static int Master() { return 0; } // MasterProc; }

   static void muteStdOstreams() {
      if (!IamMaster())
         mute_->Mute();
   }
   static void RecoverStdOstreams() { mute_->Recover(); }
   static inline double Wtime() { return MPI_Wtime(); };
   static inline void Barrier() { MPI_Barrier(Communicator_); };
   static inline void Abort(int code) { MPI_Abort(Communicator_, code); };

   static const MPI_Comm &myComm() { return Communicator_; }

   static void setLoopSplitWork(std::string f) {
      loopDivision = (f == "block") ? 1 : 0;
   }
   static const int first_element(int n) {
      // int nloc = n / size_;
      // nloc += (my_rank_==size_-1)? n%size : 0;
      // return my_rank_ * nloc;
      return (loopDivision == 0) ? my_rank_
                                 : my_rank_ * static_cast<int>(n / size_);
   }
   static const int next_element(int n) {
      return (loopDivision == 0) ? size_ : 1;
   }
   static const int last_element(int n) {
      if (loopDivision == 0)
         return n;
      else {
         return (my_rank_ == size_ - 1)
                    ? n
                    : (my_rank_ + 1) * static_cast<int>(n / size_);
      }
   }

   ~MPIcf();

   static inline long Send(double &a, int who) {
      return WSend(&a, 1, who, MPI_TAG<double>::TAG, Communicator_, rq);
   }
   static inline long Send(long &a, int who) {
      return WSend(&a, 1, who, MPI_TAG<long>::TAG, Communicator_, rq);
   }

   static inline long Recv(double &a, int who) {
      return WRecv(&a, 1, who, MPI_TAG<double>::TAG, Communicator_, rq);
   }
   static inline long Recv(long &a, int who) {
      return WRecv(&a, 1, who, MPI_TAG<long>::TAG, Communicator_, rq);
   }

   template <typename T> static inline long Bcast(T &a, int who, int size);
   template <typename T>
   static inline long Bcast(T &a, int who, int size, const MPI_Comm &comm);
   //   return WBcast(&a, 1, who,comm);}
   // static inline long Bcast(int &a, int who)   {return WBcast(&a, 1,
   // who,Communicator_);} static inline long Bcast(int &a, int who, const
   // MPI_Comm& comm)   {
   //   return WBcast(&a, 1, who,comm);}
   // static inline long Bcast(long &a, int who)   {return WBcast(&a, 1,
   // who,Communicator_);} static inline long Bcast(long &a, int who,const
   // MPI_Comm& comm )   {
   //   return WBcast(&a, 1, who,comm);}

   // template<typename T>
   // static inline long Bcast(T &a, int who); {return WBcast(&a, 1,
   // who,Communicator_);}
   template <typename T> static inline long Bcast(const KN<T> &a, int who);
   template <typename T>
   static inline long Bcast(const KN<T> &a, int who, const MPI_Comm &comm);

   template <typename T>
   static inline void Reduce(const T &, T &, int, const MPI_Op &, int);
   template <typename T>
   static inline void Reduce(const KN<T> &, KN<T> &, const MPI_Op &, int);

   template <typename T>
   static inline void AllReduce(const T &, T &, const MPI_Op &);
   template <typename T>
   static inline void AllReduce(const T &, T &, int, const MPI_Op &);
   template <typename T>
   static inline void AllReduce(const KN<T> &, KN<T> &, const MPI_Op &);
   template <typename T>
   static inline void AllReduce(const T *, T *, int, const MPI_Op &);

   template <typename T>
   static inline void Scan(const T &, T &, const MPI_Op &);

   static inline void Dup(MPI_Comm &newcomm) {
      MPI_Comm_dup(MPI_COMM_WORLD, &newcomm);
   }
   static inline void Split(MPI_Comm *newcomm, int color) {
      int key = my_rank_;
      MPI_Comm_split(MPI_COMM_WORLD, color, key, newcomm);
   }
};

#include "cfmpi.tpp"

#endif
