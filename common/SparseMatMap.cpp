#include "SparseMatMap.hpp"


// void multiply( int N, const std::map<std::pair<int,int>,double>& A, const Rn& rhs, Rn& b){
//   b = 0.;
//   assert(b.size() == rhs.size() && b.size() == N);
//   auto itA = A.begin();
//   while(itA != A.end()) {
//     b(itA->first.first) += itA->second * rhs(itA->first.second);
//     itA++;
//   }
// }

void multiply( int N, int M, const std::map<std::pair<int,int>,double>& A, const Rn& rhs, Rn& b){
  b = 0.;
  assert(b.size() == N && rhs.size() == M);
  auto itA = A.begin();
  while(itA != A.end()) {
    b(itA->first.first) += itA->second * rhs(itA->first.second);
    itA++;
  }
}

void multiply( int N, const std::map<std::pair<int,int>,double>& AA, const std::map<std::pair<int,int>,double>& BB, std::map<std::pair<int,int>,double>& C){
  C.clear();
  SparseMatrixRC<double> A(N,N,AA);
  SparseMatrixRC<double> B(N,N,BB);

  KN<int> xb(N, -1);
  KN<double> x(N, 0.);
  long ibot = 2*max(AA.size(), BB.size());

  KN<int> JC(ibot, 0);

  int ip = 0;

  for(int i=0; i<N; ++i) {    //loop over line of A

    int ic = ip;              // the line we are computing

    for(int jp=A.p[i]; jp<A.p[i+1];++jp) {   // loop over row i of A
      int j  = A.j[jp];                         // get the column

      for(int kp=B.p[j]; kp<B.p[j+1];++kp) {     // loop over line j of B (corresponding to the column of element in A)

        int k = B.j[kp];             // get column of element in row j of B

        if( xb(k) != i){            // first time we have column k element
          JC(ip) = k;                // save column element ip
          ip += 1;
          xb(k) = i;
          x(k) = A.a[jp] * B.a[kp];
        }
        else{
          x(k) += A.a[jp] * B.a[kp];
        }
      }
    }
    if(ip > ibot) assert(0);
    for(int vp= ic; vp<ip; ++vp) {
      int v = JC(vp);
      C[make_pair(i,v)] = x(v);
    }
  }
}

void multiply(const SparseMatrixRC<double>& A, const SparseMatrixRC<double>& B, std::map<std::pair<int,int>,double>& C){
  C.clear();
  int r = B.M;
  int p = A.N;

  KN<int> xb(r, -1);
  KN<double> x(r, 0.);
  long ibot = 2*max(A.nbcoef, B.nbcoef);

  KN<int> JC(ibot, 0);

  int ip = 0;

  for(int i=0; i<p; ++i) {    //loop over line of A

    int ic = ip;              // the line we are computing

    for(int jp=A.p[i]; jp<A.p[i+1];++jp) {   // loop over row i of A
      int j  = A.j[jp];                         // get the column

      for(int kp=B.p[j]; kp<B.p[j+1];++kp) {     // loop over line j of B (corresponding to the column of element in A)

        int k = B.j[kp];             // get column of element in row j of B

        if( xb(k) != i){            // first time we have column k element
          JC(ip) = k;                // save column element ip
          ip += 1;
          xb(k) = i;
          x(k) = A.a[jp] * B.a[kp];
        }
        else{
          x(k) += A.a[jp] * B.a[kp];
        }
      }
    }
    if(ip > ibot) assert(0);
    for(int vp= ic; vp<ip; ++vp) {
      int v = JC(vp);
      C[make_pair(i,v)] = x(v);
    }
  }
}









//
