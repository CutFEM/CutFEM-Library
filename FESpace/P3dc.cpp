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

  // void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
  //   for (int i = 0; i < 3; ++i) v[i] = 1;
  // }


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


void TypeOfFE_P3dcLagrange2d::FB(const What_d whatd, const Element & K,
			     const Rd & P, RNMK_ & val) const
{

  R l[]={1.-P.x-P.y,P.x,P.y};

  assert(val.N() >= E::nv+2*E::ne+E::nf);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

    if (whatd & Fop_D0)  {

      f0[0]  = l[0]*(3.0*l[0]-1.0)*(3.0*l[0]-2.0)*0.5;
      f0[1]  = l[1]*(3.0*l[1]-1.0)*(3.0*l[1]-2.0)*0.5;
      f0[2]  = l[2]*(3.0*l[2]-1.0)*(3.0*l[2]-2.0)*0.5;

      f0[3]  = 9.0*l[1]*l[2]*(3.0*l[1]-1.0)*0.5;
      f0[4]  = 9.0*l[1]*l[2]*(3.0*l[2]-1.0)*0.5;

      f0[5]  = 9.0*l[0]*l[2]*(3.0*l[2]-1.0)*0.5;
      f0[6]  = 9.0*l[0]*l[2]*(3.0*l[0]-1.0)*0.5;

      f0[7]  = 9.0*l[0]*l[1]*(3.0*l[0]-1.0)*0.5;
      f0[8]  = 9.0*l[0]*l[1]*(3.0*l[1]-1.0)*0.5;

      f0[9] = 27.0*l[0]*l[1]*l[2];
    }

    if( whatd & Fop_D1 ) {
      R2 Dl[3];
      K.Gradlambda(Dl);

      R l0[3]={ l[0],(3.0*l[0]-1.0),(3.0*l[0]-2.0)};
      R l1[3]={ l[1],(3.0*l[1]-1.0),(3.0*l[1]-2.0)};
      R l2[3]={ l[2],(3.0*l[2]-1.0),(3.0*l[2]-2.0)};
      int d=0;
      for( int der=op_dx; der<=op_dy; der++,d++) {

        RN_ f0x(val('.',0,der));

        f0x[0] = 0.5 * (Dl[0][d]*l0[1]*l0[2] + l0[0]*3.0*Dl[0][d]*l0[2] + l0[0]*l0[1]*3.0*Dl[0][d]);
        f0x[1] = 0.5 * (Dl[1][d]*l1[1]*l1[2] + l1[0]*3.0*Dl[1][d]*l1[2] + l1[0]*l1[1]*3.0*Dl[1][d]);
        f0x[2] = 0.5 * (Dl[2][d]*l2[1]*l2[2] + l2[0]*3.0*Dl[2][d]*l2[2] + l2[0]*l2[1]*3.0*Dl[2][d]);

        f0x[3] = 9.0 * (Dl[1][d]*l[2]*l1[1] + l[1]*Dl[2][d]*l1[1] + l[1]*l[2]*3.0*Dl[1][d]) * 0.5;
        f0x[4] = 9.0 * (Dl[1][d]*l[2]*l2[1] + l[1]*Dl[2][d]*l2[1] + l[1]*l[2]*3.0*Dl[2][d]) * 0.5;

        f0x[5] = 9.0 * (Dl[0][d]*l[2]*l2[1] + l[0]*Dl[2][d]*l2[1] + l[0]*l[2]*3.0*Dl[2][d]) * 0.5;
        f0x[6] = 9.0 * (Dl[0][d]*l[2]*l0[1] + l[0]*Dl[2][d]*l0[1] + l[0]*l[2]*3.0*Dl[0][d]) * 0.5;

        f0x[7] = 9.0 * (Dl[0][d]*l[1]*l0[1] + l[0]*Dl[1][d]*l0[1] + l[0]*l[1]*3.0*Dl[0][d]) * 0.5;
        f0x[8] = 9.0 * (Dl[0][d]*l[1]*l1[1] + l[0]*Dl[1][d]*l1[1] + l[0]*l[1]*3.0*Dl[1][d]) * 0.5;

        f0x[9] = 27.0 * (Dl[0][d]*l[1]*l[2] + l[0]*Dl[1][d]*l[2] + l[0]*l[1]*Dl[2][d]);
      }


      if (whatd & Fop_D2) {

      	R2 Dl0[3]={ Dl[0],(3.0*Dl[0]),(3.0*Dl[0])};
      	R2 Dl1[3]={ Dl[1],(3.0*Dl[1]),(3.0*Dl[1])};
      	R2 Dl2[3]={ Dl[2],(3.0*Dl[2]),(3.0*Dl[2])};
      	int d=0;
      	int op_dder[] = {op_dxx, op_dxy, op_dyy};

      	for( int d1=0; d1<2; d1++) {
      	  for( int d2=d1; d2<2; d2++) {

      	    RN_ f0x(val('.',0,op_dder[d]));

      	    f0x[0] = 0.5 * (
      			    Dl[0][d1]*Dl0[1][d2]*l0[2] + Dl[0][d1]*l0[1]*Dl0[2][d2] +
      			    Dl0[0][d2]*3.0*Dl[0][d1]*l0[2] + l0[0]*3.0*Dl[0][d1]*Dl0[2][d2] +
      			    Dl0[0][d2]*l0[1]*3.0*Dl[0][d1] + l0[0]*Dl0[1][d2]*3.0*Dl[0][d1]
      			    );
      	    f0x[1] = 0.5 * (
      			    Dl[1][d1]*Dl1[1][d2]*l1[2] + Dl[1][d1]*l1[1]*Dl1[2][d2] +
      			    Dl1[0][d2]*3.0*Dl[1][d1]*l1[2] + l1[0]*3.0*Dl[1][d1]*Dl1[2][d2] +
      			    Dl1[0][d2]*l1[1]*3.0*Dl[1][d1] + l1[0]*Dl1[1][d2]*3.0*Dl[1][d1]
      			    );
      	    f0x[2] = 0.5 * (
      			    Dl[2][d1]*Dl2[1][d2]*l2[2] + Dl[2][d1]*l2[1]*Dl2[2][d2] +
      			    Dl2[0][d2]*3.0*Dl[2][d1]*l2[2] + l2[0]*3.0*Dl[2][d1]*Dl2[2][d2] +
      			    Dl2[0][d2]*l2[1]*3.0*Dl[2][d1] + l2[0]*Dl2[1][d2]*3.0*Dl[2][d1]
      			    );

      	    f0x[3] = 9.0 * (
      			    Dl[1][d1]*Dl[2][d2]*l1[1] + Dl[1][d1]*l[2]*Dl1[1][d2] +
      			    Dl[1][d2]*Dl[2][d1]*l1[1] + l[1]*Dl[2][d1]*Dl1[1][d2] +
      			    Dl[1][d2]*l[2]*3.0*Dl[1][d1] + l[1]*Dl[2][d2]*3.0*Dl[1][d1]
      			    ) * 0.5;
      	    f0x[4] = 9.0 * (
      			    Dl[1][d1]*Dl[2][d2]*l2[1] + Dl[1][d1]*l[2]*Dl2[1][d2] +
      			    Dl[1][d2]*Dl[2][d1]*l2[1] + l[1]*Dl[2][d1]*Dl2[1][d2] +
      			    Dl[1][d2]*l[2]*3.0*Dl[2][d1] + l[1]*Dl[2][d2]*3.0*Dl[2][d1]
      			    ) * 0.5;


      	    f0x[5] = 9.0 * (
      			    Dl[0][d1]*Dl[2][d2]*l2[1] + Dl[0][d1]*l[2]*Dl2[1][d2] +
      			    Dl[0][d2]*Dl[2][d1]*l2[1] + l[0]*Dl[2][d1]*Dl2[1][d2] +
      			    Dl[0][d2]*l[2]*3.0*Dl[2][d1] + l[0]*Dl[2][d2]*3.0*Dl[2][d1]
      			    ) * 0.5;
      	    f0x[6] = 9.0 * (
      			    Dl[0][d1]*Dl[2][d2]*l0[1] + Dl[0][d1]*l[2]*Dl0[1][d2] +
      			    Dl[0][d2]*Dl[2][d1]*l0[1] + l[0]*Dl[2][d1]*Dl0[1][d2] +
      			    Dl[0][d2]*l[2]*3.0*Dl[0][d1] + l[0]*Dl[2][d2]*3.0*Dl[0][d1]
      			    ) * 0.5;

      	    f0x[7] = 9.0 * (
      			    Dl[0][d1]*Dl[1][d2]*l0[1] + Dl[0][d1]*l[1]*Dl0[1][d2] +
      			    Dl[0][d2]*Dl[1][d1]*l0[1] + l[0]*Dl[1][d1]*Dl0[1][d2] +
      			    Dl[0][d2]*l[1]*3.0*Dl[0][d1] + l[0]*Dl[1][d2]*3.0*Dl[0][d1]
      			    ) * 0.5;
      	    f0x[8] = 9.0 * (
      			    Dl[0][d1]*Dl[1][d2]*l1[1] + Dl[0][d1]*l[1]*Dl1[1][d2] +
      			    Dl[0][d2]*Dl[1][d1]*l1[1] + l[0]*Dl[1][d1]*Dl1[1][d2] +
      			    Dl[0][d2]*l[1]*3.0*Dl[1][d1] + l[0]*Dl[1][d2]*3.0*Dl[1][d1]
      			    ) * 0.5;

      	    f0x[9] = 27.0 * (
      			     Dl[0][d1]*Dl[1][d2]*l[2] + Dl[0][d1]*l[1]*Dl[2][d2] +
      			     Dl[0][d2]*Dl[1][d1]*l[2] + l[0]*Dl[1][d1]*Dl[2][d2] +
      			     Dl[0][d2]*l[1]*Dl[2][d1] + l[0]*Dl[1][d2]*Dl[2][d1]
      			     );

      	  d++;
      	  }
      	}

      }
    }

  }



static TypeOfFE_P3dcLagrange2d  P3dc_2d;
GTypeOfFE<Mesh2> & P3dcLagrange2d(P3dc_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P3dc=P3dc_2d;
