#include "FESpace.hpp"


// P1
class TypeOfFE_P1dcLagrange2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 1;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P1dcLagrange2d(): GTypeOfFE<Mesh2>(3, 1, Data, 3, 3, alpha_Pi_h) {

    static const R2 Pt[3] = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};

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

const int TypeOfFE_P1dcLagrange2d::nbNodeOnItem[4] = {1,0,0,0};
int TypeOfFE_P1dcLagrange2d::Data[] = {
  6, 6, 6,    // we use the face because we want discontinuous element
  0, 1, 2,    // the number of the df on  the node
  0, 0, 0,    // the node of the df
  0, 1, 2,    // which are de df on sub FE
  0, 0, 1, 0, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  3           // end_dfcomp
};

double TypeOfFE_P1dcLagrange2d::alpha_Pi_h[] = {1. ,1. ,1.};


void TypeOfFE_P1dcLagrange2d::FB(const What_d whatd, const Element & K,  const R2 & P,RNMK_ & val) const
{
  R l[] = {1.-P.sum(),P.x,P.y};

  assert(val.N() >= Element::nv);
  assert(val.M()==1 );
  double s = 1.;///pow(K.mesure(),1./4);//1./sqrt(K.mesure());
  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = s*l[0];
    f0[1] = s*l[1];
    f0[2] = s*l[2];
  }

  if (whatd & Fop_D1) {
    R2 Dl[3];
    K.Gradlambda(Dl);
    if (whatd & Fop_dx)  {
      RN_ f0x(val('.',0,op_dx));
      f0x[0] = s*Dl[0].x;
      f0x[1] = s*Dl[1].x;
      f0x[2] = s*Dl[2].x;
    }

    if (whatd & Fop_dy) {
      RN_ f0y(val('.',0,op_dy));
      f0y[0] = s*Dl[0].y;
      f0y[1] = s*Dl[1].y;
      f0y[2] = s*Dl[2].y;
    }
  }
}

// P1 3D
class TypeOfFE_P1dcLagrange3d : public GTypeOfFE<Mesh3> {

  typedef   Mesh3 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 1;
  static const int ndf = 4;
  static int Data[];
  static double alpha_Pi_h[];

  TypeOfFE_P1dcLagrange3d(): GTypeOfFE<Mesh3>(4, 1, Data, 4, 4, alpha_Pi_h) {

    static const R3 Pt[4] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.) };

    for(int i=0;i<ndf;++i) {
      Pt_Pi_h[i] = Pt[i];
      ipj_Pi_h[i] = IPJ(i,i,0);
    }
  }


  void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P1dcLagrange3d::nbNodeOnItem[4] = {1,0,0,0};

// tre trick is to use the face because we want discontinuous element
// the the node on the face will have 3 df, one for each ndfonVertex
// like P3 on edge
int TypeOfFE_P1dcLagrange3d::Data[] = {
  14, 14, 14, 14,    // we use the volume because we want discontinuous element
  0, 1, 2, 3,    // the number of the df on  the node
  0, 0, 0, 0,    // the node of the df
  0, 1, 2, 3,    // which are de df on sub FE
  0, 0, 0, 1, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  4           // end_dfcomp
};

double TypeOfFE_P1dcLagrange3d::alpha_Pi_h[] = {1. ,1. ,1., 1.};


void TypeOfFE_P1dcLagrange3d::FB(const What_d whatd, const Element & K,  const R3 & P,RNMK_ & val) const
{
  R l[]={1.-P.sum(),P.x,P.y,P.z};

  assert(val.N() >=Element::nv);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));
  if (whatd & Fop_D0) {
    f0[0] = l[0];
    f0[1] = l[1];
    f0[2] = l[2];
    f0[3] = l[3];
  }
  if (whatd & Fop_D1) {
    R3 Dl[4];
    K.Gradlambda(Dl);

    if (whatd & Fop_dx) {
      RN_ f0x(val('.',0,op_dx));
      f0x[0] = Dl[0].x;
      f0x[1] = Dl[1].x;
      f0x[2] = Dl[2].x;
      f0x[3] = Dl[3].x;
    }

    if (whatd & Fop_dy) {
      RN_ f0y(val('.',0,op_dy));
      f0y[0] = Dl[0].y;
      f0y[1] = Dl[1].y;
      f0y[2] = Dl[2].y;
      f0y[3] = Dl[3].y;
    }

    if (whatd & Fop_dz) {
      RN_ f0z(val('.',0,op_dz));
      f0z[0] = Dl[0].z;
      f0z[1] = Dl[1].z;
      f0z[2] = Dl[2].z;
      f0z[3] = Dl[3].z;
    }
  }
}


// P1
class TypeOfFE_P1dcTaylor2d : public GTypeOfFE<Mesh2> {

  typedef   Mesh2 Mesh;
  typedef  typename Mesh::Element  E;
  static const int nbNodeOnItem[4];
public:

  static const int k = 1;
  static const int ndf = (k + 2) * (k + 1) / 2;
  static int Data[];
  static double alpha_Pi_h[];
  static const double h_eps;
  TypeOfFE_P1dcTaylor2d(): GTypeOfFE<Mesh2>(3, 1, Data, 5, 5, alpha_Pi_h) {

    static const R2 Pt[5] = { R2(1./3., 1./3.)
                             ,R2(1./3.+h_eps, 1./3.), R2(1./3.-h_eps, 1./3.)
                             ,R2(1./3., 1./3.+h_eps), R2(1./3., 1./3.-h_eps)
                           };

    // u_h(x_G)
    Pt_Pi_h[0] = Pt[0];
    ipj_Pi_h[0] = IPJ(0,0,0);
    int ii=1;
    for(int i=1;i<3;++i) {
      Pt_Pi_h[ii] = Pt[ii];
      ipj_Pi_h[ii] = IPJ(i,ii,0);
      ii++;
      Pt_Pi_h[ii] = Pt[ii];
      ipj_Pi_h[ii] = IPJ(i,ii,0);
      ii++;
    }
  }


void FB(const What_d ,const Element & ,const Rd &, RNMK_ &) const;
} ;

const int TypeOfFE_P1dcTaylor2d::nbNodeOnItem[4] = {1,0,0,0};
const double TypeOfFE_P1dcTaylor2d::h_eps = 1e-10;
// tre trick is to use the face because we want discontinuous element
// the the node on the face will have 3 df, one for each ndfonVertex
// like P3 on edge
int TypeOfFE_P1dcTaylor2d::Data[] = {
  6, 6, 6,    // we use the face because we want discontinuous element
  0, 1, 2,    // the number of the df on  the node
  0, 0, 0,    // the node of the df
  0, 1, 2,    // which are de df on sub FE
  0, 0, 1, 0, // nb node on what
  0,          // for each compontant $j=0,N-1$ it give the sub FE associated
  0,          // begin_dfcomp
  3           // end_dfcomp
};

double TypeOfFE_P1dcTaylor2d::alpha_Pi_h[] = {1. ,1./(2.*h_eps), -1./(2.*h_eps), 1./(2.*h_eps), -1./(2.*h_eps)};


void TypeOfFE_P1dcTaylor2d::FB(const What_d whatd, const Element & K,  const R2 & Phat,RNMK_ & val) const
{

  R2 X = K(Phat);
  R2 G = K.centroid();
  R meas = K.mesure();

  assert(val.N() >= Element::nv);
  assert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

  if (whatd & Fop_D0) {
    f0[0] = 1.;
    f0[1] = (X.x - G.x);
    f0[2] = (X.y - G.y);
  }

  if (whatd & Fop_D1) {
    if (whatd & Fop_dx)  {
      RN_ f0x(val('.',0,op_dx));
      f0x[0] = 0.;
      f0x[1] = 1.;
      f0x[2] = 0.;
    }

    if (whatd & Fop_dy) {
      RN_ f0y(val('.',0,op_dy));
      f0y[0] = 0.;
      f0y[1] = 0.;
      f0y[2] = 1.;
    }
  }
}


static TypeOfFE_P1dcLagrange2d  P1dc_2d;
static TypeOfFE_P1dcLagrange3d  P1dc_3d;
static TypeOfFE_P1dcTaylor2d  P1dcTaylor_2d;
GTypeOfFE<Mesh2> & P1dcLagrange2d(P1dc_2d);
GTypeOfFE<Mesh3> & P1dcLagrange3d(P1dc_3d);
GTypeOfFE<Mesh2> & P1dcTaylor2d(P1dcTaylor_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P1dc=P1dc_2d;
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P1dc=P1dc_3d;
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P1dcTaylor=P1dcTaylor_2d;
