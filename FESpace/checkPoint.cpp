#include "checkPoint.hpp"


void checkPoint::save(const KN_<double>& v, string filename) {
  std::ofstream f;
  f.open(filename.c_str(), std::ofstream::out);
  f << v.size();
  for(int i=0;i<v.size();++i) f << v(i);
  f.close();
}
