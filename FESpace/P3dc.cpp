#include "FESpace.hpp"


// P1
class TypeOfFE_P3dcLagrange2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  // static const int nbNodeOnItem[4];
public:

  static const int k = 3;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];
  static const int nn[10][3];
  static const int aa[10][3];
  static const int ff[10];
  static const int il[10];
  static const int jl[10];
  static const int kl[10];
  // dof, dim Im, Data, coefInterp, nbPt, coeff
  TypeOfFE_P3dcLagrange2d(): GTypeOfFE<Mesh2>(10, 1, Data, 10, 10, alpha_Pi_h) {

    static const R2 Pt[10] = {R2(0 / 3., 0 / 3.), R2(3 / 3., 0 / 3.), R2(0 / 3., 3 / 3.),
      R2(2 / 3., 1 / 3.), R2(1 / 3., 2 / 3.),
      R2(0 / 3., 2 / 3.), R2(0 / 3., 1 / 3.),
      R2(1 / 3., 0 / 3.), R2(2 / 3., 0 / 3.),
      R2(1 / 3., 1 / 3.)};

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
      ipj_Pi_h[i] = IPJ(i,i,0);
    }
  }

  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

// const int TypeOfFE_P3dcLagrange2d::nbNodeOnItem[4] = {1,2,1,0};
int TypeOfFE_P3dcLagrange2d::Data[] = {
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6,   // the support number  of the node of the df
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,   // the number of the df on  the node
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // the node of the df
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,   // which are de df on sub FE
  0, 0, 1, 0,        // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  10           // end_dfcomp
};

double TypeOfFE_P3dcLagrange2d::alpha_Pi_h[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

// const R2 TypeOfFE_P3dcLagrange2d::G(1. / 3., 1. / 3.);
// const R TypeOfFE_P3dcLagrange2d::cshrink = 1 - 1e-2;
// const R TypeOfFE_P3dcLagrange2d::cshrink1 = 1. / TypeOfFE_P3dcLagrange2d::cshrink;


void TypeOfFE_P3dcLagrange2d::FB(const What_d whatd, const Element & K,
			     const Rd & PHat, RNMK_ & val) const
{
  R2 A(K[0]), B(K[1]), C(K[2]);
  R l0 = 1. - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
  R L[3] = {l0 * k, l1 * k, l2 * k};
  assert(val.N() >= 10);
  assert(val.M()==1 );

  int p[10] = {};

  for (int i = 0; i < 10; ++i) {
      p[i] = i;
  }
  val=0;
  RN_ f0(val('.',0,op_id));

    if (whatd & Fop_D0)  {
      RN_ f0(val('.', 0, op_id));

      for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R f = 1. / ff[df];

          for (int i = 0; i < k; ++i) {
              f *= L[nn[df][i]] - aa[df][i];
          }

          f0[pdf] = f;
      }
    }

    if( whatd & Fop_D1 || whatd & Fop_D2) {
      R2 D[] = {K.H(0) * k, K.H(1) * k, K.H(2) * k};

      if( whatd & Fop_D1){

        for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R fx = 0., fy = 0., f = 1. / ff[df];

          for (int i = 0; i < k; ++i) {
            int n = nn[df][i];
            R Ln = L[n] - aa[df][i];
            fx = fx * Ln + f * D[n].x;
            fy = fy * Ln + f * D[n].y;
            f = f * Ln;
          }
          val(pdf, 0, op_dx) = fx;
          val(pdf, 0, op_dy) = fy;

        }
      }

      if (whatd & Fop_D2) {
        for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R fx = 0., fy = 0., f = 1. / ff[df];
          R fxx = 0., fyy = 0., fxy = 0.;

          for (int i = 0; i < k; ++i) {
            int n = nn[df][i];
            R Ln = L[n] - aa[df][i];
            fxx = fxx * Ln + 2. * fx * D[n].x;
            fyy = fyy * Ln + 2. * fy * D[n].y;
            fxy = fxy * Ln + fx * D[n].y + fy * D[n].x;
            fx = fx * Ln + f * D[n].x;
            fy = fy * Ln + f * D[n].y;
            f = f * Ln;
          }
          val(pdf, 0, op_dxx) = fxx;
          val(pdf, 0, op_dyy) = fyy;
          val(pdf, 0, op_dxy) = fxy;
        }
      }
    }

  }


  const int TypeOfFE_P3dcLagrange2d::nn[10][3] = {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {1, 1, 2},
                                                {1, 2, 2}, {0, 2, 2}, {0, 0, 2}, {0, 0, 1},
                                                {0, 1, 1}, {0, 1, 2}};
  const int TypeOfFE_P3dcLagrange2d::aa[10][3] = {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 0},
                                                {0, 0, 1}, {0, 0, 1}, {0, 1, 0}, {0, 1, 0},
                                                {0, 0, 1}, {0, 0, 0}};
  const int TypeOfFE_P3dcLagrange2d::ff[10] = {6, 6, 6, 2, 2, 2, 2, 2, 2, 1};
  const int TypeOfFE_P3dcLagrange2d::il[10] = {3, 0, 0, 0, 0, 1, 2, 2, 1, 1};
  const int TypeOfFE_P3dcLagrange2d::jl[10] = {0, 3, 0, 2, 1, 0, 0, 1, 2, 1};
  const int TypeOfFE_P3dcLagrange2d::kl[10] = {0, 0, 3, 1, 2, 2, 1, 0, 0, 1};


static TypeOfFE_P3dcLagrange2d  P3dc_2d;
GTypeOfFE<Mesh2> & P3dcLagrange2d(P3dc_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P3dc=P3dc_2d;
