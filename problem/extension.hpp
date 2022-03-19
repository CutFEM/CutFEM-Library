#ifndef _EXTENSION_HPP
#define _EXTENSION_HPP


#include "../FESpace/macroElement.hpp"
#include "baseCutProblem.hpp"
#include<set>




class Extension {
public:
  const int good = 0, extension = 1, exhaust = 2;

  const QuadratureFormular1d& QF  = QF_GaussLegendre2;
  const QuadratureFormular2d& QFK = QuadratureFormular_T_5;

  set<int> dof2rm;
  std::map<std::pair<int,int>,double> S, St;
  std::map<std::pair<int,int>,double> P, Pt;

  CutFEM<Mesh2>& problem;
  std::map<const FESpace2*,std::pair<int,int>> mapIdx0_K; // space -> idxK_begin , idxK_end
  std::map<const FESpace2*,const GMacro* > mapMacro; // space -> macro_class

public:
  std::map<std::pair<int ,int>, int> element_edge_handle; //(id_element, id_edge) => what_to_do

  Extension(CutFEM<Mesh2>& pp) : problem(pp) {
    for(auto it =problem.mapIdx0_K.begin();it != problem.mapIdx0_K.end();++it){
      std:pair<int,int> v(it->second, it->second+it->first->NbElement());
      mapIdx0_K[it->first] = v;
      mapMacro[it->first] = nullptr;
    }
    element_edge_handle.clear();
  }

  public:
    void tag_extension_edges(const MacroElement& macro) {return tag_extension_edges(macro, macro.Vh);};
    void tag_extension_edges(const MacroElement& macro, const FESpace2& Vh);
    void tag_extension_edges(const MacroElement& macro, const CHyperFace& b);
    void tag_extension_edges(const MacroElement& macro, const CHyperFace& ed, const CBorder& bo);

    void tag_exhaust_edges(const MacroElement& macro) ;
    void solve(string solverName = "mumps");
    void solve_weak_extension(string solverName = "mumps");

    void do_extension();
  private:
    void do_extension_edge(const std::map<std::pair<int,int>,int>::const_iterator& it);
    void do_extension_element(const SmallElement&);

    void do_extension_P0  (const FElement2& Ks, const FElement2& Kf, int ic);
    void do_weak_extension_P0  (const FElement2& Ks, const FElement2& Kf, int ic);

  //   void do_extension_P1dc  (const FElement2& Ks, const FElement2& Kf, int ic);
    void do_extension_RT0 (const FElement2& Ks, const FElement2& Kf, int id_e, int ic);
    void do_weak_extension_RT0 (const FElement2& Ks, const FElement2& Kf, int id_e, int ic);

  //   void do_extension_BDM1(const FElement2& Ks, const FElement2& Kf, int id_e, int ic);
  //   void do_extension_RT1(const FElement2& Ks, const FElement2& Kf, int id_e, int ic);
  //
  //   void evaluate_dofP1dc(const FElement2& FKs, const FElement2& FKf, Rnm& val, int ic) ;
    void evaluate_dofRT0(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) ;
    void evaluate_dofBDM1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val);
  //   void evaluate_dofRT1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val);


  // int map_to_global(const FESpace2& vh, int k) const {
  //   auto it = mapIdx0_K.find(vh);
  //   return it->second.first + k;
  // }
  // int map_to_local(const FESpace2& vh, int k) const {
  //   auto it = mapIdx0_K.find(vh);
  //   return  k - it->second.first;
  // }
  int k_begin(const FESpace2& vh) const {
    auto it = mapIdx0_K.find(&vh);
    assert(it != mapIdx0_K.end());
    return it->second.first;
  }
  const FESpace2& get_space_from_idxK(int k, int& i0) const {
    for(auto it = mapIdx0_K.begin(); it!=mapIdx0_K.end();++it){
      i0 = it->second.first;
      if(k>= it->second.first && k < it->second.second) return *(it->first);
    }
    assert(0);
  }
  const GMacro& get_macro(const FESpace2& vh) const {
    auto it = mapMacro.find(&vh);
    assert(it != mapMacro.end());
    return *(it->second);
  }

private:
public:
    void make_S();
    void make_weak_extension();
    void precond(Rn& rhs);
    void removeDF( int N, std::map<std::pair<int,int>,double>& A, Rn& b);
    void reconstruct(Rn& b);
public:
    void erase_rows_to_fix_RT0();
    void erase_rows_to_fix2_RT0();
    void erase_rows_to_fix_average_RT0();

    void erase_rows_to_fix_BDM1();
// friend void save(const MacroElement & macro, const Extension& extension) ;
};





#endif
