#ifndef _MACRO_ELEMENT_HPP
#define _MACRO_ELEMENT_HPP


#include "CutFESpace.hpp"
#include<set>

class Extension;

struct MElement {
public:

  int idx_root_element;
  vector<int> idx_element;
  vector<std::pair<int,int>> inner_edge; // idx_Element, idx_edge

  // MElement(int idx_root=-1, int idx_K=-1, int idx_edge=-1){
  //   idx_root_element = idx_root;
  //   idx_element.push_back(idx_root);
  //   int kk =(idx_root_element < idx_K)? idx_root_element : idx_K;
  //   this->add(idx_K, std::make_pair(kk, idx_edge));
  // }
  MElement(int idx_root = -1){
    idx_root_element = idx_root;
    idx_element.push_back(idx_root);
  }

  void add(int idx_K, std::pair<int,int> ke) {
    idx_element.push_back(idx_K);
    // int kk =(idx_root_element < idx_K)? idx_root_element : idx_K;
    inner_edge.push_back(ke);
  }

  void print() const {
    std::cout << " Root idx   : \t " << idx_root_element << std::endl;
    std::cout << " small element : \n";
    for(auto it = inner_edge.begin(); it != inner_edge.end();++it) {
      std::cout << it->first << " with edge " << it->second << std::endl;
    }
  }
};

struct SmallElement {
  int index;
  int index_root;
  int chain_position;
  int idx_edge_to_root;

  SmallElement(int idx = -1, int idx_root = -1) :index(idx), index_root(idx_root), chain_position(0), idx_edge_to_root(-1) {}
  void setRoot(int i) { index_root = i;}
  void setChainPosition(int i) { chain_position = i;}
  void setEdgeDirection(int i) { idx_edge_to_root = i;}
};


class GMacro {
public :

const int small=-1 ;
const int good = 0, extension = 1, exhaust = 2;


map<int, MElement>     macro_element;   // idx_root -> idx_macroElement
map<int, SmallElement> small_element;  // idx_element -> idx_small_element
double tol;

GMacro() {}

int getIndexRootElement(int k) const {
  auto it = small_element.find(k);
  if(it == small_element.end()){
    return k;
  }else {
    return it->second.index_root;
  }
}
const SmallElement& getSmallElement(int k) const {
  auto it = small_element.find(k);
  assert(it != small_element.end());
  return it->second;
}
bool isRootFat(int k) const {
  return (small_element.find(k) == small_element.end()) && (k!=-1);
}

};

class MacroElementSurface : public GMacro {
  typedef typename Interface2::FaceIdx Face;
public:
  const Interface2& interface;

  MacroElementSurface(const Interface2& , const double) ;
  void find_small_element() ;
  void find_root_element()  ;
  int check_direction(const int, const int, int&);
};

class MacroElement : public GMacro {

public:
  
  const FESpace2& Vh;

  R tol;
  int nb_element_0, nb_element_1;

  MacroElement(const FESpace2& vh, const double C);



private:
  void find_small_element();
  void find_root_element();


  friend class Extension;
};




#endif
