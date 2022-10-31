#include "FESpace.hpp"


// P1 Polynomial basis {1, x}
class TypeOfFE_P0Polynomial1d : public  GTypeOfFE<Mesh1> {

  typedef   Mesh1 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 0;
  static const int ndf = (k + 1);
  static int Data[];
  static double alpha_Pi_h[];

  // dof, dim Im, Data, coefInterp, nbPt, coeff
  TypeOfFE_P0Polynomial1d():GTypeOfFE<Mesh1>(1, 1, Data, 1, 1, alpha_Pi_h){

    GTypeOfFE<Mesh1>::basisFctType = BasisFctType::P0poly;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R1 Pt[1] = {R1(0.5)};

    for(int i=0;i<ndf;++i) Pt_Pi_h[i] = Pt[i];
    for(int i=0;i<1;  ++i) ipj_Pi_h[i] = IPJ(i,i,0);

  }



  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P0Polynomial1d::nbNodeOnItem[4] = {1,0,0,0};
int TypeOfFE_P0Polynomial1d::Data[] = {
  2,    // the support number  of the node of the df
  0,    // the number of the df on  the node
  0,    // the node of the df
  0,    // which are de df on sub FE
  0, 1, 0, 0, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  1           // end_dfcomp
};
//
double TypeOfFE_P0Polynomial1d::alpha_Pi_h[] = {1.};
void TypeOfFE_P0Polynomial1d::FB(const What_d whatd, const Element & K,
			       const R1 & P,RNMK_ & val) const
{
  assert(val.N() >= 1);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = 1.;
  }

  if (whatd & Fop_D1)
    {
      if (whatd & Fop_dx)
  	{
  	  RN_ f0x(val('.',0,op_dx));
  	  f0x[0] = 0.;
  	}
  }
}
static TypeOfFE_P0Polynomial1d  P0_Polynomial_1d;
GTypeOfFE<Mesh1> & P0Polynomial1d(P0_Polynomial_1d);
template<> GTypeOfFE<Mesh1> & DataFE<Mesh1>::P0Poly=P0_Polynomial_1d;


// P0 2D
class TypeOfFE_P0Lagrange2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 0;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P0Lagrange2d(): GTypeOfFE<Mesh2>(1, 1, Data, 1, 1, alpha_Pi_h) {

    GTypeOfFE<Mesh2>::basisFctType = BasisFctType::P0;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R2 Pt[1] = {R2(1./3, 1./3)};

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

const int TypeOfFE_P0Lagrange2d::nbNodeOnItem[4] = {0,0,1,0};
int TypeOfFE_P0Lagrange2d::Data[] = {
  6,    // the support number  of the node of the df
  0,    // the number of the df on  the node
  0,    // the node of the df
  0,    // which are de df on sub FE
  0, 0, 1, 0, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  1           // end_dfcomp
};
double TypeOfFE_P0Lagrange2d::alpha_Pi_h[] = {1.};


void TypeOfFE_P0Lagrange2d::FB(const What_d whatd, const Element & K,
			       const R2 & P,RNMK_ & val) const
{
  assert(val.N() >= 1);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = 1.;
  }

}
static TypeOfFE_P0Lagrange2d  P0_2d;
GTypeOfFE<Mesh2> & P0Lagrange2d(P0_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P0=P0_2d;



// P0 SCALE 2D
class TypeOfFE_P0ScLagrange2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 0;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  // static double alpha_Pi_h[];

  TypeOfFE_P0ScLagrange2d(): GTypeOfFE<Mesh2>(1, 1, Data, 1, 1) {

    GTypeOfFE<Mesh2>::basisFctType = BasisFctType::P0;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R2 Pt[1] = {R2(1./3, 1./3)};

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
      ipj_Pi_h[i] = IPJ(i,i,0);
    }
  }


  void get_Coef_Pi_h(const GbaseFElement<Mesh> & K, KN_<double> &v) const;


  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P0ScLagrange2d::nbNodeOnItem[4] = {0,0,1,0};
int TypeOfFE_P0ScLagrange2d::Data[] = {
  6,    // the support number  of the node of the df
  0,    // the number of the df on  the node
  0,    // the node of the df
  0,    // which are de df on sub FE
  0, 0, 1, 0, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  1           // end_dfcomp
};
// double TypeOfFE_P0ScLagrange2d::alpha_Pi_h[] = {1.};
void TypeOfFE_P0ScLagrange2d::get_Coef_Pi_h(const GbaseFElement<Mesh> & K, KN_<double> &v) const {
  const Element &T = K.T;
  v[0] = T.measure();
}

void TypeOfFE_P0ScLagrange2d::FB(const What_d whatd, const Element & K,
			       const R2 & P,RNMK_ & val) const
{
  assert(val.N() >= 1);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = 1./K.measure();
  }

}
static TypeOfFE_P0ScLagrange2d  P0sc_2d;
GTypeOfFE<Mesh2> & P0ScLagrange2d(P0sc_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P0sc=P0sc_2d;


// P0 3D
class TypeOfFE_P0Lagrange3d : public GTypeOfFE<Mesh3> {

  typedef   Mesh3 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 0;
  static const int ndf = 1;//(k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P0Lagrange3d(): GTypeOfFE<Mesh3>(1, 1, Data, 1, 1, alpha_Pi_h) {
    GTypeOfFE<Mesh3>::basisFctType = BasisFctType::P0;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R3 Pt[1] = {R3(1./4, 1./4, 1./4)};

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

const int TypeOfFE_P0Lagrange3d::nbNodeOnItem[4] = {0,0,0,1};
int TypeOfFE_P0Lagrange3d::Data[] = {
  14,   // the support number  of the node of the df
  0,    // the number of the df on  the node
  0,    // the node of the df
  0,    // which are de df on sub FE
  0, 0, 0, 1, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  1           // end_dfcomp
};
double TypeOfFE_P0Lagrange3d::alpha_Pi_h[] = {1.};


void TypeOfFE_P0Lagrange3d::FB(const What_d whatd, const Element & K,
			       const R3 & P,RNMK_ & val) const
{
  assert(val.N() >= 1);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = 1.;
  }

}


static TypeOfFE_P0Lagrange3d  P0_3d;
GTypeOfFE<Mesh3> & P0Lagrange3d(P0_3d);
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P0=P0_3d;


// P0 3D HEXA
class TypeOfFE_P0LagrangeQ3d : public GTypeOfFE<MeshHexa> {

  typedef   MeshHexa Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 0;
  static const int ndf = 1;//(k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P0LagrangeQ3d(): GTypeOfFE<MeshHexa>(1, 1, Data, 1, 1, alpha_Pi_h) {
    GTypeOfFE<Mesh>::basisFctType = BasisFctType::P0;
    GTypeOfFE<Mesh>::polynomialOrder = k;

    static const R3 Pt[1] = {R3(1./2, 1./2, 1./2)};

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

const int TypeOfFE_P0LagrangeQ3d::nbNodeOnItem[4] = {0,0,0,1};
int TypeOfFE_P0LagrangeQ3d::Data[] = {
  14,   // the support number  of the node of the df
  0,    // the number of the df on  the node
  0,    // the node of the df
  0,    // which are de df on sub FE
  0, 0, 0, 1, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  1           // end_dfcomp
};
double TypeOfFE_P0LagrangeQ3d::alpha_Pi_h[] = {1.};


void TypeOfFE_P0LagrangeQ3d::FB(const What_d whatd, const Element & K,
			       const R3 & P,RNMK_ & val) const
{
  assert(val.N() >= 1);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = 1.;
  }

}


static TypeOfFE_P0LagrangeQ3d  P0Q_3d;
GTypeOfFE<MeshHexa> & P0LagrangeQ3d(P0Q_3d);
template<> GTypeOfFE<MeshHexa> & DataFE<MeshHexa>::P0=P0Q_3d;
