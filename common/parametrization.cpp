#include "parametrization.hpp"


LinearParametrization::LinearParametrization(const std::vector<double> &x, const std::vector<double> & y) {
  this->init(x, y);
}


void LinearParametrization::init(const std::vector<double> &x, const std::vector<double> & y) {

  int n = x.size()-1;
  nb_interval = n;

  std::vector<double> b(n);
  for(int i = 1; i < n; ++i) b[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);

  myLineSet = new LineSet[n];
  for(int i = 0; i < n; ++i) {
    mySplineSet[i].a = y[i];
    mySplineSet[i].b = b[i];
    mySplineSet[i].x = x[i];
  }
  return;
}
std::vector
