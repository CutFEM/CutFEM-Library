#include "FESpace.hpp"



// P1 Polynomial basis {1, x}
class TypeOfFE_P2Polynomial1d : public  GTypeOfFE<Mesh1> {

  typedef   Mesh1 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 2;
  static const int ndf = (k + 1);
  static int Data[];
  static double alpha_Pi_h[];

  // dof, dim Im, Data, coefInterp, nbPt, coeff
  TypeOfFE_P2Polynomial1d(): GTypeOfFE<Mesh1>(3, 1, Data, 7, 3, alpha_Pi_h) {
    GTypeOfFE<Mesh>::basisFctType = BasisFctType::P2poly;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R1 Pt[3] = {R1(0.), R1(1.), R1(1./2)};

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
    }

    ipj_Pi_h[0] = IPJ(0,0,0);
    for(int i=0;i<3;++i) {
      ipj_Pi_h[i+1] = IPJ(1,i,0);
    }
    for(int i=0;i<3;++i) {
      ipj_Pi_h[i+4] = IPJ(2,i,0);
    }
  }



  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P2Polynomial1d::nbNodeOnItem[4] = {1,1,0,0};
int TypeOfFE_P2Polynomial1d::Data[] = {
  0, 1, 2,    // the support number  of the node of the df
  0, 0, 0,    // the number of the df on  the node
  0, 1, 2,    // the node of the df
  0, 1, 2,    // which are de df on sub FE
  1, 1, 0, 0, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  3           // end_dfcomp
};

double TypeOfFE_P2Polynomial1d::alpha_Pi_h[] = {1., -3., 4., -1., 2., -4., 2. };

void TypeOfFE_P2Polynomial1d::FB(const What_d whatd, const Element & K, const R1 & P,RNMK_ & val) const {
  assert(K[0].x < K[1].x);
  assert(0<=P.x && P.x <=1);
  R l[]={1,P.x, P.x*P.x};

  assert(val.N() >=Element::nv);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0)
    {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
    }
  if (whatd & Fop_D1)
    {
      if (whatd & Fop_dx)
  	{
  	  RN_ f0x(val('.',0,op_dx));
  	  f0x[0] = 0;
  	  f0x[1] = 1./K.mesure();
      f0x[2] = 2*l[1]*f0x[1];
  	}
  }
}


static TypeOfFE_P2Polynomial1d  P2_Polynomial_1d;
GTypeOfFE<Mesh1> & P2Polynomial1d(P2_Polynomial_1d);
template<> GTypeOfFE<Mesh1> & DataFE<Mesh1>::P2Poly=P2_Polynomial_1d;


// P2 - 2D
class TypeOfFE_P2Lagrange2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 2;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P2Lagrange2d(): GTypeOfFE<Mesh2>(6, 1, Data, 6, 6, alpha_Pi_h) {
    GTypeOfFE<Mesh>::basisFctType = BasisFctType::P2;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R2 Pt[6] = {R2(0., 0.), R2(1., 0.), R2(0., 1.),
			     R2(1./2, 1./2), R2(0.,1./2), R2(1./2,0)};

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
      ipj_Pi_h[i] = IPJ(i,i,0);
    }
  }


  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P2Lagrange2d::nbNodeOnItem[4] = {1,1,0,0};
int TypeOfFE_P2Lagrange2d::Data[] = {
  0, 1, 2, 3, 4, 5,   // the support number  of the node of the df
  0, 0, 0, 0, 0, 0,   // the number of the df on  the node
  0, 1, 2, 3, 4, 5,   // the node of the df
  0, 1, 2, 3, 4, 5,   // which are de df on sub FE
  1, 1, 0, 0,        // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  6           // end_dfcomp
};

double TypeOfFE_P2Lagrange2d::alpha_Pi_h[] = {1. ,1. ,1., 1., 1., 1.};


void TypeOfFE_P2Lagrange2d::FB(const What_d whatd, const Element & K,
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


static TypeOfFE_P2Lagrange2d  P2_2d;
GTypeOfFE<Mesh2> & P2Lagrange2d(P2_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P2=P2_2d;







// P2 - 3D
class TypeOfFE_P2Lagrange3d : public GTypeOfFE<Mesh3> {

  typedef   Mesh3 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 2;
  static const int ndf = 10;//(k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P2Lagrange3d(): GTypeOfFE<Mesh3>(10, 1, Data, 10, 10, alpha_Pi_h) {
    GTypeOfFE<Mesh>::basisFctType = BasisFctType::P2;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R3 Pt[10] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.),
			      R3(1./2, 0., 0.), R3(0., 1./2, 0.), R3(0., 0., 1./2),
			      R3(1./2, 1./2, 0.), R3(1./2, 0.,1./2), R3(0., 1./2, 1./2)};

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
      ipj_Pi_h[i] = IPJ(i,i,0);
    }
  }


  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P2Lagrange3d::nbNodeOnItem[4] = {1,1,0,0};
int TypeOfFE_P2Lagrange3d::Data[] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  // the support number  of the node of the df
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // the number of the df on  the node
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  // the node of the df
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  // which are de df on sub FE
  1, 1, 0, 0,        // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  10           // end_dfcomp
};

double TypeOfFE_P2Lagrange3d::alpha_Pi_h[] = {1. ,1. ,1., 1., 1., 1., 1., 1., 1., 1.};


void TypeOfFE_P2Lagrange3d::FB(const What_d whatd, const Element & K,
			       const R3 & P,RNMK_ & val) const
{
  R l[]={1.-P.sum(),P.x,P.y,P.z};

  assert(val.N() >=E::nv+E::ne);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));
  if (whatd & Fop_D0) {
    int k=0;
    for(int i=0;i<E::nv;++i)
      f0[k++] = l[i]*(2*l[i]-1.);
    for(int i=0;i<E::ne;++i)
      f0[k++] = 4.*l[E::nvedge[i][0]]*l[E::nvedge[i][1]];
  }

  if (whatd & (Fop_D1|Fop_D2)) {

    R3 Dl[4];
    R l4[4]={ (4*l[0]-1),(4*l[1]-1),(4*l[2]-1),(4*l[3]-1)};
    K.Gradlambda(Dl);

    if( whatd & (Fop_D1) ) {
      RN_ f0x(val('.',0,op_dx));
      RN_ f0y(val('.',0,op_dy));
      RN_ f0z(val('.',0,op_dz));

      int k=0;
      for(int i=0;i<E::nv;++i,++k) {
	f0x[k] = Dl[i].x*l4[i];
	f0y[k] = Dl[i].y*l4[i];
	f0z[k] = Dl[i].z*l4[i];
      }

      for(int i=0;i<E::ne;++i,++k) {
	int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	f0x[k] = 4*(Dl[i1].x*l[i0] + Dl[i0].x*l[i1]) ;
	f0y[k] = 4*(Dl[i1].y*l[i0] + Dl[i0].y*l[i1]) ;
	f0z[k] = 4*(Dl[i1].z*l[i0] + Dl[i0].z*l[i1]) ;
      }

      assert(k==10);
    }

    if (whatd & Fop_D2) {
      RN_ f0xx(val('.',0,op_dxx));
      RN_ f0yy(val('.',0,op_dyy));
      RN_ f0zz(val('.',0,op_dzz));
      RN_ f0xy(val('.',0,op_dxy));
      RN_ f0xz(val('.',0,op_dxz));
      RN_ f0yz(val('.',0,op_dyz));


      int k=0;
      for(int i=0;i<E::nv;++i,++k)
	{
	  f0xx[k] = 4.*Dl[i].x*Dl[i].x;
	  f0yy[k] = 4.*Dl[i].y*Dl[i].y;
	  f0zz[k] = 4.*Dl[i].z*Dl[i].z;
	  f0xy[k] = 4.*Dl[i].x*Dl[i].y;
	  f0xz[k] = 4.*Dl[i].x*Dl[i].z;
	  f0yz[k] = 4.*Dl[i].y*Dl[i].z;
	}
      for(int i=0;i<E::ne;++i,++k)
	{
	  int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	  f0xx[k] = 8.*Dl[i0].x*Dl[i1].x;
	  f0yy[k] = 8.*Dl[i0].y*Dl[i1].y;
	  f0zz[k] = 8.*Dl[i0].z*Dl[i1].z;
	  f0xy[k] = 4.*(Dl[i0].x*Dl[i1].y+ Dl[i1].x*Dl[i0].y);
	  f0xz[k] = 4.*(Dl[i0].x*Dl[i1].z+ Dl[i1].x*Dl[i0].z);
	  f0yz[k] = 4.*(Dl[i0].y*Dl[i1].z+ Dl[i1].y*Dl[i0].z);
	}
      assert(k==10);
    }
  }
}


static TypeOfFE_P2Lagrange3d  P2_3d;
GTypeOfFE<Mesh3> & P2Lagrange3d(P2_3d);
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P2=P2_3d;
