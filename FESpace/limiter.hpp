
#include "FESpace.hpp"
#include "macroElement.hpp"




namespace limiter {

  // Only for discontinuous Lagrange
  template<typename Mesh>
  void extendToMacroP0(const FunFEM<Mesh>& uh, Rn& u_new, const MacroElement<Mesh>& macro) {
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename Mesh::Rd Rd;

    const GFESpace<Mesh>& Wh(*uh.Vh);
    for (auto& p : macro.macro_element) {

      // compute the value of the dof of the elements in the
      // macro element.
      // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K \cap \Omega|}
      const MElement& MK(p.second);
      int n_element = MK.size();

      // loop over the elements that has to be changed
      for(auto s : MK.idx_element ) {
        const FElement& FK(Wh[s]);
        assert(FK.NbDoF() == FK.NbNode());
        double area_total = MK.area_total_;

        // loop over the dof (in that case just node evaluation)
        for(int df=0; df<FK.NbNode(); ++df) {
          Rd P = FK.Pt(df);
          int df_glb = FK.loc2glb(df);
          double val = 0;

          for(auto k : MK.idx_element ) {
            double s = macro.get_area(k);
            val += s * uh.eval(k, P);
          }
          val = val / area_total;
          u_new(df_glb) = val;
        }
      }

    }
  }

  // Only for discontinuous Lagrange
  template<typename Mesh>
  void extendToMacroP1(const FunFEM<Mesh>& uh, Rn& u_new, const MacroElement<Mesh>& macro) {
    typedef typename Mesh::Element Element;

    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename Mesh::Rd Rd;

    const QF& qf(*QF_Simplex<typename FElement::RdHat>(2));
    const GFESpace<Mesh>& Wh(*uh.Vh);
    FunFEM<Mesh> fun_u_new(Wh, u_new);

    for (auto& p : macro.macro_element) {

      // compute the value of the dof of the elements in the
      // macro element.
      // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K \cap \Omega|}
      const MElement& MK(p.second);
      int n_element = MK.size();
      double area_total = MK.area_total_;

      // loop over the elements that has to be changed
      for(auto s : MK.idx_element ) {
        const FElement& FK(Wh[s]);

        // loop over the dof (in that case just node evaluation)
        for(int df=0; df<FK.NbDoF() ; ++df) {
          Rd P = FK.Pt(df);
          int df_glb = FK.loc2glb(df);
          double val = 0;

          for(auto k : MK.idx_element ) {
            double s = macro.get_area(k);
            val += s * uh.eval(k, P);
          }
          val = val / area_total;
          u_new(df_glb) = val;
        }
      }

      double C0 = 0;
      for(auto k : MK.idx_element ) {
        const FElement& FK(Wh[k]);
        const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));

        for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
          double meas = cutK.measure(it);

          for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip = cutK.mapToPhysicalElement(it, ip);
            double Cint = meas * ip.getWeight();
            C0 += Cint*(uh.eval(k, mip) - fun_u_new.eval(k,mip));
          }
        }
      }
      C0 = C0 / area_total;

      for(auto k : MK.idx_element ) {
        const FElement& FK(Wh[k]);
        for(int df=0; df<FK.NbDoF() ; ++df) {
          int df_glb = FK.loc2glb(df);
          u_new(df_glb) += C0;
        }
      }
    }
  }


};

// class Limiter {
//
//   typedef typename Mesh2::Element Element;
//   typedef typename FESpace2::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::QFB QFB;
//
//   const QF&  qf ;
//   const QFB& qfb;
//
// public:
//   Rn indicator;
//
//   Limiter() : qf(*QF_Simplex<R2>(3)), qfb(*QF_Simplex<R1>(3)){}
//
//   void KXRCF_indicator(const Fun2_h& uh, const Fun2_h& flux);
//
//   void limiting(const Fun2_h& uh, Rn& lh);
//   void min_and_max_average_neighbor(const FElement& FK, const Fun2_h&uh, double&m, double&M);
//   double average_on_Kj(const FElement& FK, const Fun2_h&uh);
//   double compute_alpha(const FElement& FK, const Fun2_h&uh, double Uj, double m, double M);
// };
//
//
// class Box_Side {
// public:
//   // a side is the part that satisfies ax+by+c < 0
//   double a, b, c;
//   Box_Side(double aa, double bb, double cc) : a(aa), b(bb), c(cc) {}
//
//   double f(const R2 A) const {
//     return (a*A.x+b*A.y+c);
//   }
//   bool is_negative(const R2 P) const {
//     return (f(P) <=0);
//   }
//   bool is_positive(const R2 P) const {
//     return (f(P) > 0);
//   }
//   bool cut_edge(const R2 A, const R2 B) const {
//     return (f(A)*f(B)< 0);
//   }
//   R2 get_cut_node(const R2 A, const R2 B) const {
//     double t = -f(A)/(f(B)-f(A));
//     return (1-t) * A + t * B;
//   }
//
// };
//
// class Filter_Box {
//   typedef typename Mesh2::Element Element;
// public:
//   static const int nv = 3;                         // triangle
//   typedef SortArray<Ubyte, nv> TriaIdx;           // the vertices of a triangle of the cut:
//
//   std::vector<Box_Side> sides_;
//   std::vector<R2>      vertices_;
//   std::vector<TriaIdx> triangles_;               // idx of the vertices
//
//   R2 minXY, maxXY;
//
//   Filter_Box(const R2 mminXY, const R2 mmaxXY) : minXY(mminXY), maxXY(mmaxXY){
//     sides_.push_back(Box_Side( 0., -1.,  minXY.y));
//     sides_.push_back(Box_Side( 1.,  0., -maxXY.x));
//     sides_.push_back(Box_Side( 0.,  1., -maxXY.y));
//     sides_.push_back(Box_Side(-1.,  0.,  minXY.x));
//   }
//
//   bool contain(const R2 P) const{
//     for(auto it=sides_.begin(); it != sides_.end();++it) {
//       if(it->is_positive(P)) return false;
//     }
//     return true;
//   }
//
//   bool intersect_with(const typename Mesh2::Element& K) const{
//     // if we find one node inside the box it is enough to know
//     for(int i=0; i<3;++i) {
//       if(this->contain((R2) K[i])) return true;
//     }
//     return false;
//   }
//
//   void cut_triangle(const typename Mesh2::Element& K) ;
//
//   void print(std::string filename = "newTriangle.dat") const ;
//
// };
//
// class Filter {
// public:
//
//
// };
