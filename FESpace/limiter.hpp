
#include "FESpace.hpp"
#include "CutFESpace.hpp"
class Limiter {

  typedef typename Mesh2::Element Element;
  typedef typename FESpace2::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::QFB QFB;

  const QF&  qf ;
  const QFB& qfb;

public:
  Rn indicator;

  Limiter() : qf(*QF_Simplex<R2>(3)), qfb(*QF_Simplex<R1>(3)){}

  void KXRCF_indicator(const Fun2_h& uh, const Fun2_h& flux);

  void limiting(const Fun2_h& uh, Rn& lh);
  void min_and_max_average_neighbor(const FElement& FK, const Fun2_h&uh, double&m, double&M);
  double average_on_Kj(const FElement& FK, const Fun2_h&uh);
  double compute_alpha(const FElement& FK, const Fun2_h&uh, double Uj, double m, double M);
};


class Box_Side {
public:
  // a side is the part that satisfies ax+by+c < 0
  double a, b, c;
  Box_Side(double aa, double bb, double cc) : a(aa), b(bb), c(cc) {}

  double f(const R2 A) const {
    return (a*A.x+b*A.y+c);
  }
  bool is_negative(const R2 P) const {
    return (f(P) <=0);
  }
  bool is_positive(const R2 P) const {
    return (f(P) > 0);
  }
  bool cut_edge(const R2 A, const R2 B) const {
    return (f(A)*f(B)< 0);
  }
  R2 get_cut_node(const R2 A, const R2 B) const {
    double t = -f(A)/(f(B)-f(A));
    return (1-t) * A + t * B;
  }

};

class Filter_Box {
  typedef typename Mesh2::Element Element;
public:
  static const int nv = 3;                         // triangle
  typedef SortArray<Ubyte, nv> TriaIdx;           // the vertices of a triangle of the cut:

  std::vector<Box_Side> sides_;
  std::vector<R2>      vertices_;
  std::vector<TriaIdx> triangles_;               // idx of the vertices

  R2 minXY, maxXY;

  Filter_Box(const R2 mminXY, const R2 mmaxXY) : minXY(mminXY), maxXY(mmaxXY){
    sides_.push_back(Box_Side( 0., -1.,  minXY.y));
    sides_.push_back(Box_Side( 1.,  0., -maxXY.x));
    sides_.push_back(Box_Side( 0.,  1., -maxXY.y));
    sides_.push_back(Box_Side(-1.,  0.,  minXY.x));
  }

  bool contain(const R2 P) const{
    for(auto it=sides_.begin(); it != sides_.end();++it) {
      if(it->is_positive(P)) return false;
    }
    return true;
  }

  bool intersect_with(const typename Mesh2::Element& K) const{
    // if we find one node inside the box it is enough to know
    for(int i=0; i<3;++i) {
      if(this->contain((R2) K[i])) return true;
    }
    return false;
  }

  void cut_triangle(const typename Mesh2::Element& K) ;

  void print(std::string filename = "newTriangle.dat") const ;

};

class Filter {
public:


};
