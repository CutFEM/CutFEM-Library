#include "FESpace.hpp"


// P1
class TypeOfFE_P2dcLagrange2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 2;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  // dof, dim Im, Data, coefInterp, nbPt, coeff
  TypeOfFE_P2dcLagrange2d(): GTypeOfFE<Mesh2>(6, 1, Data, 6, 6, alpha_Pi_h) {

    static const R2 Pt[6] = {R2(0., 0.), R2(1., 0.), R2(0., 1.),
			     R2(1./2, 1./2), R2(0.,1./2), R2(1./2,0)};

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
      ipj_Pi_h[i] = IPJ(i,i,0);
    }
  }

  // void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
  //   for (int i = 0; i < 3; ++i) v[i] = 1;
  // }


  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P2dcLagrange2d::nbNodeOnItem[4] = {1,1,0,0};
int TypeOfFE_P2dcLagrange2d::Data[] = {
  6, 6, 6, 6, 6, 6,   // the support number  of the node of the df
  0, 1, 2, 3, 4, 5,   // the number of the df on  the node
  0, 0, 0, 0, 0, 0,   // the node of the df
  0, 1, 2, 3, 4, 5,   // which are de df on sub FE
  0, 0, 1, 0,        // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  6           // end_dfcomp
};

double TypeOfFE_P2dcLagrange2d::alpha_Pi_h[] = {1. ,1. ,1., 1., 1., 1.};


void TypeOfFE_P2dcLagrange2d::FB(const What_d whatd, const Element & K,
				const R2 & P,RNMK_ & val) const
{
  R l[]={1.-P.sum(),P.x,P.y};

  assert(val.N() >=E::nv+E::ne);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0)  {
    int k=0;
    for(int i=0;i<E::nv;++i)
      f0[k++] = l[i]*(2*l[i]-1.);
    for(int i=0;i<E::ne;++i)
      f0[k++] = 4.*l[E::nvedge[i][0]]*l[E::nvedge[i][1]];
  }

  if (whatd & (Fop_D1|Fop_D2)) {
    R2 Dl[3];
    R l4[3]={ (4*l[0]-1),(4*l[1]-1),(4*l[2]-1)};

    K.Gradlambda(Dl);
    RN_ f0x(val('.',0,op_dx));
    RN_ f0y(val('.',0,op_dy));
    int k=0;
    for(int i=0;i<E::nv;++i,++k) {
      f0x[k] = Dl[i].x*l4[i];
      f0y[k] = Dl[i].y*l4[i];
    }
    for(int i=0;i<E::ne;++i,++k) {
      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
      f0x[k] = 4*(Dl[i1].x*l[i0] + Dl[i0].x*l[i1]) ;
      f0y[k] = 4*(Dl[i1].y*l[i0] + Dl[i0].y*l[i1]) ;
    }
    assert(k==6);

    if (whatd & Fop_D2) {
      RN_ f0xx(val('.',0,op_dxx));
      RN_ f0yy(val('.',0,op_dyy));
      RN_ f0xy(val('.',0,op_dxy));

      k=0;
      for(int i=0;i<E::nv;++i,++k) {
    	f0xx[k] = 4.*Dl[i].x*Dl[i].x;
    	f0yy[k] = 4.*Dl[i].y*Dl[i].y;
    	f0xy[k] = 4.*Dl[i].x*Dl[i].y;
      }
      for(int i=0;i<E::ne;++i,++k) {
    	int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
    	f0xx[k] = 8.*Dl[i0].x*Dl[i1].x;
    	f0yy[k] = 8.*Dl[i0].y*Dl[i1].y;
    	f0xy[k] = 4.*(Dl[i0].x*Dl[i1].y+ Dl[i1].x*Dl[i0].y);
      }
    }
  }
}


static TypeOfFE_P2dcLagrange2d  P2dc_2d;
GTypeOfFE<Mesh2> & P2dcLagrange2d(P2dc_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P2dc=P2dc_2d;
