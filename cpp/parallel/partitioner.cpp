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
#include "partitioner.hpp"
#include "../problem/solver.hpp"

void MetisMatrixOrdering::buildGraphOfMatrix(
    const std::map<std::pair<int, int>, double> &M) {

   // beginAdjacencyNodes_ = new idx_t[nbOfVertices_+1];
   // adjacencyArray_ = new idx_t[M.size()];

   beginAdjacencyNodes_.resize(nbOfVertices_ + 1);
   adjacencyArray_.resize(M.size());

   int k                   = 0;
   beginAdjacencyNodes_[0] = 0;

   for (typename std::map<std::pair<int, int>, double>::const_iterator q =
            M.begin();
        q != M.end(); ++q) {
      int i = q->first.first - iStart_;
      int j = q->first.second;
      if (i + iStart_ == j)
         continue;
      beginAdjacencyNodes_[i + 1] = k + 1;
      adjacencyArray_[k]          = q->first.second;
      ++k;
   }
}

void MetisMatrixOrdering::performOrdering() {

   // perm_ = new idx_t[nbOfVertices_];
   // iperm_ = new idx_t[nbOfVertices_];
   perm_.resize(nbOfVertices_);
   iperm_.resize(nbOfVertices_);

   int msg = METIS_NodeND(&nbOfVertices_, beginAdjacencyNodes_, adjacencyArray_,
                          nullptr, nullptr, perm_, iperm_);
}

void ParMetisMatrixOrdering::performOrdering() {

   idx_t numflag = 0;
   idx_t wgtflag = 0;
   idx_t option[3];
   option[0]    = 0;
   idx_t ncon   = 1;
   idx_t nparts = MPIcf::size();
   idx_t dbglvl = 5;
   perm_.resize(indexDistribution_[MPIcf::size()]);
   perm_ = 0.0;
   iperm_.resize(indexDistribution_[MPIcf::size()]);
   iperm_ = 0.0;

   int *pointerToIperm_ = perm_ + indexDistribution_[MPIcf::my_rank()];

   idx_t *size   = new idx_t[2 * MPIcf::size()];
   idx_t edgecut = 0;
   KN<real_t> tpwgts;
   tpwgts = 1.0;
   KN<real_t> ubvec;
   ubvec = 1.05;

   MPI_Comm theCommunicator;
   MPIcf::Dup(theCommunicator);

   int msg = ParMETIS_V32_NodeND(
       indexDistribution_, beginAdjacencyNodes_, adjacencyArray_, nullptr,
       &numflag, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &dbglvl,
       pointerToIperm_, size, &theCommunicator);

   MPIcf::AllReduce(perm_, iperm_, MPI_SUM);

   // int msg = ParMETIS_V3_PartKway(indexDistribution_, beginAdjacencyNodes_,
   // adjacencyArray_, 			       nullptr, nullptr, &wgtflag, 			       &numflag, &ncon, &nparts,
   // tpwgts, ubvec, option, 			       & edgecut, perm_,  &theCommunicator);

   delete[] size;
}

void ParMetisMatrixOrdering::computeIndexBeginEnd() {

   MPIcf::Scan(nbOfVertices_, iEnd_, MPI_SUM);
   iStart_ = iEnd_ - nbOfVertices_;
}

void ParMetisMatrixOrdering::setIndexDistribution() {

   KN<int> idxDist(MPIcf::size() + 1);
   idxDist = 0.0;
   KN<int> idxDistRecv(MPIcf::size() + 1);
   idxDistRecv = 0.0;
   for (int i = MPIcf::my_rank() + 1; i <= MPIcf::size(); ++i)
      idxDist(i) = nbOfVertices_;

   MPIcf::AllReduce(idxDist, idxDistRecv, MPI_SUM);

   // indexDistribution_ = new idx_t[idxDist.size()];
   indexDistribution_.resize(idxDist.size());

   for (int i = 0; i <= MPIcf::size(); ++i)
      indexDistribution_[i] = idxDistRecv[i];
}

void ConsecutiveMatrixOrdering::performOrdering() {

   // perm_ = new idx_t[nbOfVertices_];

   // int k=0, i0=-1;
   // for (typename std::map<std::pair<int,int>,double>::const_iterator
   // q=M.begin();
   //      q != M.end(); ++q)
   //   {
   //     int i =q->first.first;

   //     if(i0 != i) {
   // 	perm_[i] = k;
   // 	++k;
   //     }

   //   }
   // assert(k==M.size());
}
