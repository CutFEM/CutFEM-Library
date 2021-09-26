#ifndef _MACRO_ELEMENT_HPP
#define _MACRO_ELEMENT_HPP


#include "CutFESpace.hpp"

class MElement {
public:

  int idx_root_element;
  vector<int> idx_element;
  vector<std::pair<int,int>> inner_edge;

  // MElement(){};
  MElement(int idx_root, int idx_K, int idx_edge){
    idx_root_element = idx_root;
    idx_element.push_back(idx_root);
    this->add(idx_K, idx_edge);
  }

  void add(int idx_K, int idx_edge) {
    idx_element.push_back(idx_K);
    inner_edge.push_back(std::make_pair(idx_K, idx_edge));
  }

  void print() const {
    std::cout << " Root idx   : \t " << idx_root_element << std::endl;
    std::cout << " small element : \n";
    for(auto it = inner_edge.begin(); it != inner_edge.end();++it) {
      std::cout << it->first << " with edge " << it->second << std::endl;
    }
  }
};

class MacroElement {

const QuadratureFormular1d& QF = QF_GaussLegendre2;

public:

  enum{small=-1};

  vector<int> small_or_fat_K; // all the element [-1 or idx]
  vector<int> idx_small_K_temp;   // list of just the small element
  vector<std::pair<int,int>> idx_small_K;   // list of just the small element
  vector<std::pair<int,int>> edges_to_stabilize;
  map<int,int> chain_position;
  set<int> dof2rm;

  // maybe only need one map <root, MElement>
  vector<MElement> macro_element;
  map<int, int> idxRoot2idxMElement;

  std::map<std::pair<int,int>,double> S, St;
  std::map<std::pair<int,int>,double> P, Pt;

  const FESpace2& Vh;

  R tol;
  int nb_element_0, nb_element_1;

  MacroElement(const FESpace2& vh, const double C) : Vh(vh), small_or_fat_K(vh.nbElement){
    double h = Vh[0].T.lenEdge(0);
    double meas = Vh[0].T.mesure();
    nb_element_0 = 0;
    nb_element_1 = 0;
    tol = C * meas;

    for(int i=0;i<small_or_fat_K.size();++i) small_or_fat_K[i] = i;

    // std::cout << small_or_fat_K.size()<< std::endl;
    std::cout << " tolerance \t" << tol << std::endl;

    find_small_K();
    std::cout << " Found " << idx_small_K_temp.size() << " small elements" << std::endl;
    std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
    std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
    find_face_to_stabilize() ;
    std::cout << " Number of edge found \t " << edges_to_stabilize.size() << std::endl;

    // for(int k=0; k<idx_small_K.size(); ++k){
    //
    //   int idx = idx_small_K[k].first;
    //   std::cout << idx << "\t" <<   Vh.whichDomain(idx) << "\t"
    //             << chain_position[idx] << std::endl;
    //
    // }


  }

  bool isFat(int k) const {
    return (small_or_fat_K[k] != small);
  }

  bool isRootFat(int k) const {
    return (small_or_fat_K[k] == k);
  }

  void find_small_K() {
    for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
      if(!Vh.isCut(k)) continue;

      const FElement2& FK(Vh[k]);
      const int kb = Vh.Th(FK.T);
      const int domain = FK.whichDomain();

      CutData2 cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
      const Partition2& cutK =  Partition2(FK.T, cutData);  // build the cut
      ElementSignEnum the_part = cutK.what_part(domain);
      double areaCut = cutK.mesure(domain);

      if(areaCut < tol) {
        small_or_fat_K[k] = small;   // -1 to small elements
        idx_small_K_temp.push_back(k);
        chain_position[k] = 0;
        if(domain == 0) { nb_element_0++;}else{ nb_element_1++;}
      }
    }
  }

  void find_face_to_stabilize() {
    int ii = 0;
    while (idx_small_K_temp.size() > 0) {

      int nb_of_small_K_left = idx_small_K_temp.size();
      for (int i=nb_of_small_K_left-1;i>=0;--i) {

        int k = idx_small_K_temp[i];
        const FElement2& FK(Vh[k]);
        int k_back = Vh.Th(FK.T);
        int the_domain = FK.whichDomain();

        int kn, ie;
        bool found_fat_neigh = false;
        for(int ifac = 0; ifac < 3; ++ifac) {    //loop over the edges / faces

          int ifacn = ifac;
          int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
          if(kn_back == -1) continue;
          int kn_tmp = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
          if(kn_tmp ==-1) continue;


          if(isRootFat(kn_tmp)) {
            kn = kn_tmp;
            ie = ifac;
            chain_position[k] = 1;
            found_fat_neigh = true;

            idxRoot2idxMElement[kn] = macro_element.size();
            macro_element.push_back(MElement(kn, k, ie));
            break;
          }
          else if(isFat(kn_tmp)) {
            if(chain_position[k] == 0 ||  chain_position[k] > chain_position[kn_tmp]){
              kn = kn_tmp;
              ie = ifac;
              chain_position[k] = chain_position[kn_tmp]+1;
              found_fat_neigh = true;
            }
          }

        }

        // if(isFat(kn)) {
        if(found_fat_neigh) {
          edges_to_stabilize.push_back(std::make_pair(k, ie));
          idx_small_K.push_back(std::make_pair(k, small_or_fat_K[kn]));
          small_or_fat_K[k] = small_or_fat_K[kn];
          idx_small_K_temp.erase(idx_small_K_temp.begin()+i); // remove the element.


          if(chain_position[k] != 1){
            int idx_root = small_or_fat_K[kn];
            int idx_ME = idxRoot2idxMElement[idx_root];
            macro_element[idx_ME].add(k, ie);
          }
        }

      }
      ii++;
    }
    small_or_fat_K.clear();
  }

  void print(std::string filename = "EdgeToStabilize.dat") {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);

    for(std::vector<std::pair<int,int>>::iterator it = edges_to_stabilize.begin(); it != edges_to_stabilize.end(); ++it){

      int k = it->first;
      int ifac = it->second;

      R2 p1 = Vh[k].T[Mesh2::Element::nvedge[ifac][0]];
      R2 p2 = Vh[k].T[Mesh2::Element::nvedge[ifac][1]];

      plot << p1 << "\n" << p2 << "\n \n";
    }
    plot.close();
  }

  void make_S(){
    typedef typename Mesh2::Partition Partition;
    typedef CutData2 CutData;

    for(int i=0; i<Vh.NbDoF(); ++i){
      S[std::make_pair(i ,i)] = 1;
      St[std::make_pair(i ,i)] = 1;
    }

    dof2rm.clear();

    for(auto it=idx_small_K.begin(); it !=idx_small_K.end();++it) {

      int idx_Ks = it->first;
      int idx_Kf = it->second;

      const FElement2& Ks(Vh[idx_Ks]);
      const FElement2& Kf(Vh[idx_Kf]);

      int ks_back = Vh.Th(Ks.T);
      int kf_back = Vh.Th(Kf.T);
      int the_domain = Ks.whichDomain();


      double val[3];
      int k = idx_Ks;
      int pos_k = chain_position[k];

      for(int e=0;e<3;++e) {     // loop over edges

        int je = e;
        int kn_back = Vh.Th.ElementAdj(ks_back, je);
        if(kn_back == -1) continue;
        int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
        if(isRootFat(kn)) continue; // good dof, common with a Fat root element

        int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
        if (pos_kn < pos_k) continue;
        if (pos_kn == pos_k && kn < k) continue;

        double C = (kn == -1)? 1. : 0.5;

        evaluate_dof(Ks, e, Kf, val);
        for(int i = Kf.dfcbegin(0); i < Kf.dfcend(0); ++i) {
          S[std::make_pair(Ks.loc2glb(e),Ks.loc2glb(e))] = 0.;        // 1st component small element
          St[std::make_pair(Ks.loc2glb(e),Ks.loc2glb(e))] = 0.;        // 1st component small element
          S[std::make_pair(Ks.loc2glb(e),Kf.loc2glb(i))] = val[i];        // 1st component small element
          St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(e))] = val[i];        // 1st component small element
        }

        dof2rm.insert(Ks.loc2glb(e));
      }
      int idx0 = 3;//Kf.dfcbegin(2);   // the pressure index
      if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) continue;
      S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
      St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
      S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
      St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;

      dof2rm.insert(Ks.loc2glb(idx0));

    }
    std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
  }

  void evaluate_dof(const FElement2& FKs, int e, const FElement2& FKf, double* val) {

    for(int i=0;i<3;++i) val[i] = 0.;
    KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
    What_d Fop = Fwhatd(0);
    double meas = FKs.T.lenEdge(e);
    // R2 normal = FKs.T.N(e);
    R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
    for(int iq=0;iq<QF.getNbrOfQuads();++iq) {

      QuadraturePoint1d ip_1d(QF[iq]);
      R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
      R2 mip_Ks   = FKs.map(ip_KsHat);
      R2 ip_KfHat = FKf.T.toKref(mip_Ks);

      FKf.BF(Fop_D0, ip_KfHat, bf);

      for(int i=0;i<3;++i) {
        val[i] += meas*ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
      }
    }
  }

  void evaluate_dof2(const FElement2& FKs, int e, const FElement2& FKf, double* val) {

    for(int i=0;i<3;++i) val[i] = 0.;
    KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
    What_d Fop = Fwhatd(0);

    for(int i=0;i<3;++i) {

      double meas = FKf.T.lenEdge(i);
      R2 normal = FKf.T.EdgeOrientation(i)*FKf.T.N(i);
      for(int iq=0;iq<QF.getNbrOfQuads();++iq) {

        QuadraturePoint1d ip_1d(QF[iq]);
        R2 ip_KfHat = FKf.T.toKref(ip_1d, i);
        R2 mip_Kf   = FKf.map(ip_KfHat);
        R2 ip_KsHat = FKs.T.toKref(mip_Kf);

        FKs.BF(Fop_D0, ip_KsHat, bf);

        val[i] += meas*ip_1d.getWeight()*(bf(e,0,0)*normal.x + bf(e,1,0)*normal.y) ;
      }
    }
  }

  void precond(std::map<std::pair<int,int>,double>& A, Rn & rhs) {
    int t0 = MPIcf::Wtime();
    int N = Vh.NbDoF();
    std::map<std::pair<int,int>,double> R;
    multiply(N, St, A, R);
    A.clear();
    multiply(N, R, S, A);
    Rn b(N);
    multiply(N, N, St, rhs, b);
    rhs = b;

    std::cout << " initial size \t" << rhs.size() << std::endl;

    // need to remove the bad dof
    removeDF(N, A, rhs);

    std::cout << " size after \t" << rhs.size() << std::endl;

    std::cout << " time precond \t" << MPIcf::Wtime() - t0 << std::endl;
  }




  void removeDF( int N, std::map<std::pair<int,int>,double>& A, Rn& b){

    std::map<std::pair<int,int>,double>  C;

    int i0=0, j=0, i=0;
    for( auto it = dof2rm.begin();it != dof2rm.end();++it) {
      int iend = *it;
      while(j < iend) {
        P [make_pair(i,j)] = 1;
        Pt[make_pair(j,i)] = 1;
        ++i, ++j;
      }
      j += 1;
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



};



























#endif
