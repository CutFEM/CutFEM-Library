#ifndef PARTITIONER_HPP_
#define PARTITIONER_HPP_


#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cassert>
#include "../util/cputime.h"
#include "metis.h"
#include "parmetis.h"
// #include "cfmpi.hpp"


// enum PartitionerType {                          // type of partitioner
//   metis = 1,
//   scotch = 2
// };


// struct GraphPart {                              // struct with input for partitioner
  

  
  
// };


// class Partitioner {                            // Abstract partitioner class
  


// };




class MatrixOrdering {
protected :
  typedef KN<int>  IndexArray;

  IndexArray perm_, iperm_;
  std::map<int, int> ipermMap_;
  
  idx_t nbOfVertices_;
  int iStart_ = 0, iEnd_ = 0;
  double time_;
  
  MatrixOrdering() :  nbOfVertices_(0),time_(0) {}
  MatrixOrdering(const int nn) : nbOfVertices_(nn),time_(0) {}
  virtual void performOrdering() = 0;


public :
  int nbOfVertices() const { return nbOfVertices_;}
  virtual int perm (const int i) const { return  perm_[i];}
  virtual int iperm(const int i) const  { return iperm_[i];}
  
  virtual void reorder(const KN<double>& rhs, KN<double>& rhsMapped) const  {
    rhsMapped.resize(rhs.size());
    for(int i=0; i<rhs.size();++i) rhsMapped(iperm_[i]) = rhs(i);
  }
  
  virtual void inverseMapp(const KN<double>& rhs, KN<double>& rhsMapped) const  {
    rhsMapped.resize(rhs.size());
    for(int i=0; i<rhs.size();++i) rhsMapped(i) = rhs(iperm_[i]);
  }
  
  int iStart() const {return iStart_;}
  int iEnd() const {return iEnd_;}
};



class ConsecutiveMatrixOrdering : public MatrixOrdering {

  void performOrdering();
public:
  ConsecutiveMatrixOrdering(const int nn, const std::map<std::pair<int,int>,double> & M) :
    MatrixOrdering(nn) {
    performOrdering();
  }
  
};




class NoOrdering : public MatrixOrdering {

public :
  NoOrdering() : MatrixOrdering() {}

  void reorder(const KN<double>& rhs, KN<double>& rhsMapped) const  {
    rhsMapped.resize(rhs.size());
    rhsMapped = rhs;
  }
  void inverseMapp(const KN<double>& rhs, KN<double>& rhsMapped) const  {
    rhsMapped.resize(rhs.size());
    rhsMapped = rhs;
  }

  
  void performOrdering() {};
  
  int perm (const int i) const { return  i;}
  int iperm(const int i) const  { return i;}


  
};


// class MetisMatrixOrdering : public MatrixOrdering {

// protected:
//   IndexArray beginAdjacencyNodes_;
//   IndexArray adjacencyArray_;
  
//   virtual void performOrdering();
//   void buildGraphOfMatrix(const std::map<std::pair<int,int>,double> & M);

// public:

//   MetisMatrixOrdering(const int nn) :
//     MatrixOrdering(nn)  {}

//   MetisMatrixOrdering(const int nn, const std::map<std::pair<int,int>,double> & M) :
//     MatrixOrdering(nn) {
//     R t0 = MPIcf::Wtime();
    
//     buildGraphOfMatrix(M);
//     performOrdering();

//     time_ = MPIcf::Wtime() - t0;
//     std::cout << " Time perfomrming reordering \t" << time_ << std::endl;
//   }
  
// };




// class ParMetisMatrixOrdering : public MetisMatrixOrdering {
  
//   IndexArray indexDistribution_;

//   void computeIndexBeginEnd();
//   void setIndexDistribution();
//   void performOrdering() ;

// public:

//     ParMetisMatrixOrdering(const int nn, const std::map<std::pair<int,int>,double> & M) :
//       MetisMatrixOrdering(nn){

//     R t0 = MPIcf::Wtime();

//     computeIndexBeginEnd();
//     setIndexDistribution();
//     buildGraphOfMatrix(M);
//     performOrdering();

//     time_ = MPIcf::Wtime() - t0;
//     std::cout << " Time perfomrming reordering \t" << time_ << std::endl;
    
//   }
// };








// class MetisPartitioner {

// protected :
//   double time_;
  
// public :
//   typedef idx_t*  IndexArray;

//   idx_t ne;                            // number of element
//   idx_t nn;                            // number of nodes in the mesh
//   IndexArray eptr;                     // starting value of nodes of element in eind
//   IndexArray eind;                     // node of element
//   IndexArray vwgt;                     // weight of element
//   IndexArray vsize;                    // size of element (not used)
//   idx_t ncommon;                       // number of nodes for connected elements
//   idx_t nparts;                        // nb of partitions
//   real_t * tpwgts;                     // weight for each partition
//   IndexArray options;
//   idx_t objval;
//   IndexArray epart;                    // partition of the element
//   IndexArray npart;                    // partition of the nodes

//   template<class Mesh>
//   MetisPartitioner(const Mesh & Th, const int np)
//     : eptr(0), eind(0), vwgt(0), vsize(0), tpwgts(0), options(0),
//       epart(0), npart(0) { buildGraph(Th,np); };
  
//   template<class Mesh>
//   void Mesh2File(const Mesh & Th);

// private:
//   template<class Mesh>
//   void buildGraph(const Mesh & Th, const int np);
// public:
//   ~MetisPartitioner() {
//     if(eptr) delete[] eptr;
//     if(eind) delete[] eind;
//     if(vwgt) delete[] vwgt;
//     if(vsize) delete[] vsize;
//     if(tpwgts) delete[] tpwgts;
//     if(options) delete[] options;
//     if(epart) delete[] epart;
//     if(npart) delete[] npart;
//   }
// };



// template<class Mesh>
// void MetisPartitioner::Mesh2File(const Mesh & Th) {

//   const int n = Th.nt;
//   const int nbdf = Th[0].nv;
//   std::ofstream meshFile;
//   meshFile.open("Th.mesh",  std::ofstream::out);
//   meshFile << n << "\t 1" << std::endl;
  
//   for(int k=0;k<Th.nt;++k) {
//     const typename Mesh::Element & K(Th[k]);
//     for(int i=0;i<nbdf;++i)
//       meshFile << Th(K[i]) << "\t";
//     meshFile << std::endl;
//   }
//   meshFile.close();

// }

// template<class Mesh>
// void MetisPartitioner::buildGraph(const Mesh & Th, const int np){
//   if( np == 1) return;
//   // std::cout << " BuildGraph with Metis " << std::endl;
//   const double cpubegin = CPUtime();

//   ne = Th.nt;
//   const int nbdf = Th[0].nv;
//   nn = Th.nv;
//   eptr = new idx_t[ne+1];
//   eind = new idx_t[ne * nbdf];
//   ncommon = Mesh::Rd::d;
//   nparts = np;
//   epart = new idx_t[ne];
//   npart = new idx_t[nn];
  
//   // const double cpuend = CPUtime();
//   // std::cout << " CPU Time For BuildGraph \t " << cpuend - cpubegin << std::endl;

//   int l=0;
//   for(int k=0;k<ne;++k) {
//     const typename Mesh::Element & K(Th[k]);
//     eptr[k] = l;
//     for(int i=0;i<nbdf;++i,++l)
//       eind[l] = Th(K[i]);
//   }
//   eptr[ne] = l;
  
//   int msg = METIS_PartMeshDual(&ne, &nn,eptr,eind,vwgt,vsize,&ncommon,&nparts,tpwgts,
//   			       options,&objval,epart,npart);

//   // std::cout << " objval \t" << objval << std::endl;
//   // std::cout << " Metis output value \t" << msg << std::endl;

// }

#endif
