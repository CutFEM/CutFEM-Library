
template<typename T>
inline
long MPIcf::Bcast(T &a, int who, int size) {
  return WBcast(&a, size, who,Communicator_);}

template<typename T>
inline
long MPIcf::Bcast(T &a, int who, int size, const MPI_Comm& comm){
  return WBcast(&a, size, who,comm);
}


template<typename T>
inline long MPIcf::Bcast(const KN<T> &a, int who)  {
  assert(&a);
  CheckContigueKN(a);
  return WBcast((T *) a, a.N(), who,Communicator_);
}
template<typename T>
inline long MPIcf::Bcast(const KN<T> &a, int who, const MPI_Comm& comm )  {
  assert(&a);
  CheckContigueKN(a);
  return WBcast((T *) a, a.N(), who, comm);
}



template <typename T>
  inline void MPIcf::Reduce(const T& myData, T& globalData, int size,
			    const MPI_Op& op, int root)
  { MPI_Reduce(&myData, &globalData, size,
	       MPI_TYPE<T>::TYPE(), op, root, Communicator_); }

template <typename T>
inline void MPIcf::Reduce(const KN<T>& myData, KN<T>& globalData,
			  const MPI_Op& op, int root)
{ MPI_Reduce((T*)(myData), (T*)globalData, myData.size(),
	     MPI_TYPE<T>::TYPE(), op, root, Communicator_); }



template <typename T>
inline void MPIcf::AllReduce(const T& myData, T& globalData,
			     const MPI_Op& op)
  { MPI_Allreduce(&myData, &globalData, 1,
		  MPI_TYPE<T>::TYPE(), op, Communicator_); }


template <typename T>
inline void MPIcf::AllReduce(const T& myData, T& globalData, int size,
			     const MPI_Op& op)
  { MPI_Allreduce(&myData, &globalData, size,
		  MPI_TYPE<T>::TYPE(), op, Communicator_); }

template <typename T>
inline void MPIcf::AllReduce(const T* myData, T* globalData, int size,
  const MPI_Op& op)
  { MPI_Allreduce(myData, globalData, size,
    MPI_TYPE<T>::TYPE(), op, Communicator_); }

template <typename T>
inline void MPIcf::AllReduce(const KN<T>& myData, KN<T>& globalData,
			     const MPI_Op& op)
  {
    CheckContigueKN(myData);
    MPI_Allreduce((T*)(myData), globalData, myData.size(),
		  MPI_TYPE<T>::TYPE(), op, Communicator_);
  }


template <typename T>
inline void MPIcf::Scan(const T& myData, T& globalData, const MPI_Op& op)
{
  MPI_Scan(&myData, &globalData, 1,
	   MPI_TYPE<T>::TYPE(), op, Communicator_);
}
