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

void eraseAndSetRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, int dof2rm, int dof2set, double val){

  std::map<std::pair<int,int>,double> C;
  // std::map<std::pair<int,int>,double> D;
  std::map<std::pair<int,int>,double> P;
  for(int i=0;i<N;++i){
    P[make_pair(i,i)] = 1;
  }

  // for( auto & p : dof2rm) {
    // int i0 = dof2rm;//p.first;
    P [make_pair(dof2set,dof2set)] = 0;
    b(dof2set) = val;//p.second;

  // }
  SparseMatrixRC<double> AA (N,N,A);
  SparseMatrixRC<double> PP (N,N,P);
  multiply(PP,AA, C);
  // SparseMatrixRC<double> CC (N,N,C);
  // multiply(CC,PP, A);
  A = C;

  // for( auto & p : dof2rm) {
    // int i0 = dof2set;//p.first;
    A [make_pair(dof2set,dof2set)] = 1;
  // }

}

void eraseAndSetRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, std::map<int, double>& dof2rm){

  // std::map<std::pair<int,int>,double> C;
  std::map<std::pair<int,int>,double> P;
  for(int i=0;i<N;++i){
    P[make_pair(i,i)] = 1;
  }

  for( auto & p : dof2rm) {
    int i0 = p.first;
    P [make_pair(i0,i0)] = 0;
    b(i0) = p.second;
  }

  SparseMatrixRC<double> AA (N,N,A);
  SparseMatrixRC<double> PP (N,N,P);
  multiply(PP,AA, A);
  // A = C;

  for( auto & p : dof2rm) {
    int i0 = p.first;
    A [make_pair(i0,i0)] = 1;
  }

}

void eraseRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, std::set<int>& dof2rm){

  std::map<std::pair<int,int>,double> C;
  std::map<std::pair<int,int>,double> P;
  std::map<std::pair<int,int>,double> Pt;

  int i = 0, j=0;
  for(auto p:dof2rm) {
    while(j<p) {
      P[make_pair(i,j)] = 1;
      Pt[make_pair(j,i)] = 1;
      j++;
      i++;
    }
    j++;
  }
  while(j<N){
    P[make_pair(i,j)] = 1;
    Pt[make_pair(j,i)] = 1;
    ++i;
    ++j;
  }

  int ndf = dof2rm.size();
  int nline = N - ndf;
  int ncol  = N;
  SparseMatrixRC<double> AA (N    ,N   ,A );
  SparseMatrixRC<double> PP (nline,ncol,P );
  SparseMatrixRC<double> PPt(ncol,nline,Pt);

  multiply(PP, AA, C);
  SparseMatrixRC<double> CC(nline,ncol,C);
  multiply(CC, PPt, A);

  Rn x(nline, 0.);
  multiply(nline, ncol, P, b, x);
  b.resize(nline);
  b = x;

}






//
