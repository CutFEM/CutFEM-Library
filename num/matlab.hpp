#ifndef MATLAB_HPP
#define MATLAB_HPP

#include <map>
#include "../common/RNM.hpp"
typedef double R;

namespace matlab {
  template<class Vector, class A>
  void Export(const Vector & a, const A& mapp,std::string filename);
  template<class Vector>
  void Export(const Vector & a, std::string filename);

  template<class A>
  void Export(std::map<std::pair<A,A>,R> & dF, std::string filename);

  template<class A>
  void Export(std::map<std::pair<int,int>,R> & dF, const A& mapp, std::string filename);


  template<class Vector>
  void loadVector(Vector & v, int nv,std::string filename);
  // template<class K>
  // void ExportM(const KNM<K> & a, std::string filename);
}


namespace matlab {
  // template<class K>
  // void matlab::ExportM(const KNM<K> & a, std::string filename) {
  //   std::ofstream plot;
  //   plot.open(filename.c_str(), std::ofstream::out);
  //   plot << std::setprecision(16);
  //   for(int i=0;i<a.N();++i) {
  //     for(int j=0;j<a.M();++j) {
  //       plot << a(i,j) << "  ";
  //     }
  //     plot << std::endl;
  //   }
  //   plot.close();
  // }

  template<class Vector>
  void loadVector(Vector & v, int nv, std::string filename) {
    std::ifstream f(filename);
    if(!f) cout << " Connot open file" <<  endl;
    cout << " Load File  \"" <<filename<<"\" in a vector " << nv <<  endl;

    v.resize(nv);
    for (int i=0;i<nv;i++) f >> v[i];
    cout << " DONE "<<  endl;
    f.close();
  }

  template<class Vector>
  void Export(const Vector & a, std::string filename) {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    plot << std::setprecision(16);
    for(int i=0;i<a.size();++i) plot << a[i] << std::endl;
    plot.close();
  }

  template<class Vector, class A>
  void Export(const Vector & a, const A& mapp,std::string filename) {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    plot << std::setprecision(16);
    for(int i=0;i<a.size();++i) plot << a(mapp.iperm(i)) << std::endl;
    plot.close();
  }

  template<class A>
  void Export(std::map<std::pair<A,A>,R> & dF, std::string filename) {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    plot << std::setprecision(16);
    for (std::map<std::pair<int,int>,R>::const_iterator q=dF.begin();
  	 q != dF.end(); ++q)
      {
  	plot << q->first.first  + 1 << "\t" << q->first.second + 1
  	     << "\t" << q->second << std::endl;
      }
    plot.close();
  }

  template<class A>
  void Export(std::map<std::pair<int,int>,R> & dF, const A& mapp, std::string filename) {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    plot << std::setprecision(18);
    for (std::map<std::pair<int,int>,R>::const_iterator q=dF.begin();
	 q != dF.end(); ++q)
      {
	plot << mapp.iperm(q->first.first)  + 1 << "\t" << mapp.iperm(q->first.second) + 1
	     << "\t" << q->second << std::endl;
      }
    plot.close();
  }





template<class Mesh>
void Export(const Mesh& mesh) {
    std::ofstream plot;
    plot.open("nodes.dat", std::ofstream::out);
    for(int i=0; i<mesh.nbVertices();++i) {
      plot << mesh(i).x << "\t" << mesh(i).y << std::endl;
    }
    plot.close();

    plot.open("elements.dat", std::ofstream::out);
    for(int k=0;k<mesh.nbElmts();++k) {
      const typename Mesh::Element & T(mesh[k]);
      for(int i =0; i<3;++i)
	plot <<  mesh(T[i]) << "\t";
      plot << "0" << std::endl;
    }
    plot.close();

    plot.open("edges.dat", std::ofstream::out);
    for(int k=0;k<mesh.nbBrdElmts();++k) {
      const typename Mesh::BorderElement & T(mesh.be(k));
      for(int i =0; i<2;++i)
	plot <<  mesh(T[i]) << "\t";
      for(int i =0; i<5;++i)
	plot <<  "0 \t";

      plot << std::endl;
    }
    plot.close();
}


}
#endif
