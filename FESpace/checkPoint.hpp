#ifndef _CHECKPOINT_HPP
#define _CHECKPOINT_HPP

#include "expression.hpp"

namespace checkPoint {
  void save(const KN_<double>& v, string filename);

  template<class Mesh>
  void load(FunFEM<Mesh>& v, string path) {
    int n;
    ifstream f(path.c_str());
    if(!f) {cerr << "Load a file to KN<double> " << path << endl; exit(1);}
    cout << " Read On file \"" << path <<"\""<<  endl;
    f >> n;
    assert(n == v.Vh->NbDoF());
    for(int i=0;i<v.size();++i) f >> v(i);
  }


};




#endif
