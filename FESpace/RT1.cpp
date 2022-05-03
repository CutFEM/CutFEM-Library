#include "FESpace.hpp"

struct InitTypeOfRTk_2d {
  int k;       // order poly on edge
  int ndfi;    // nb of internal dof
  int npe;     // nb point on edge
  int ndf;     // nb dof

  KN< R > X;    // point on edge
  // KN<R> Pi_h_coef; // 1
  KN< int > Data;    // data of TypeOfFE
  const QuadratureFormular1d QFE;
  const GQuadratureFormular<R2> &QFK;

  Vertex2 verticesHat[3];
  Triangle2 TriangleHat;

  InitTypeOfRTk_2d(int KK)
    : k(KK), ndfi((k + 1) * (k)), npe(k + 1), ndf(3 * npe + ndfi), Data(4 * ndf + 7),
      QFE(-1 + 2 * npe, npe, GaussLegendre(npe), true), QFK(QuadratureFormular_T_5) {
    // int j = 0;
    int ndfe = ndf - ndfi;    //
    int o[5];

    o[0] = 0;

    for (int i = 1; i < 6; ++i) {
      o[i] = o[i - 1] + ndf;
    }

    for (int dof = 0; dof < ndf; ++dof) {
      if (dof < ndfe) {
        int e = dof / npe;
        int n = dof % npe;
        Data[o[0] + dof] = 3 + e;
        Data[o[1] + dof] = n;
        Data[o[2] + dof] = e;
        Data[o[3] + dof] = dof;
      } else {
        int n = dof - ndfe;
        Data[o[0] + dof] = 6;
        Data[o[1] + dof] = n;
        Data[o[2] + dof] = 3;
        Data[o[3] + dof] = dof;
      }
    }

    Data[o[4] + 0] = 0;
    Data[o[4] + 1] = 1;
    Data[o[4] + 2] = 1;
    Data[o[4] + 3] = 0;
    Data[o[4] + 4] = 0;
    Data[o[4] + 5] = 0;
    Data[o[4] + 6] = ndf;    // end_dfcomp

// std::cout << Data << std::endl;
// getchar();
    (R2&) verticesHat[0] = R2::KHat[0];
    (R2&) verticesHat[1] = R2::KHat[1];
    (R2&) verticesHat[2] = R2::KHat[2];
    int iv[3] = {0,1,2};
    TriangleHat.set(verticesHat, iv, 0);

    // std::cout << TriangleHat[0] << "\n"
    // << TriangleHat[1] << "\n"
    // << TriangleHat[2] << "\n";
    //
    // std::cout << QFE << std::endl;
    // std::cout << QFK << std::endl;
  }


};


class TypeOfFE_RT1_2d : public InitTypeOfRTk_2d, public GTypeOfFE<Mesh2> {
  typedef  Mesh2 Mesh; // Define 2D mesh as Mesh
  typedef  typename Mesh::Element  E; // Define mesh element as E

 public:
  static double Pi_h_coef[];
  // bool Ortho;

  TypeOfFE_RT1_2d(): InitTypeOfRTk_2d(1), GTypeOfFE<Mesh2>(
    ndf,
    2,
    Data,
    2 * 2 * 3 * QFE.n + QFK.n * 4,    // nb coef mat interpole
    3 * QFE.n + QFK.n                 // nb P interpolation
  ){
    int kkk = 0, i = 0;
    for(int e = 0; e < 3; ++e) {
      for (int p = 0; p < QFE.n; ++p) {
        R2 A(TriangleHat[Element::nvedge[e][0]]);
        R2 B(TriangleHat[Element::nvedge[e][1]]);

        ipj_Pi_h[kkk++] = IPJ(2 * e, i, 0);
        ipj_Pi_h[kkk++] = IPJ(2 * e, i, 1);
        ipj_Pi_h[kkk++] = IPJ(2 * e + 1, i, 0);
        ipj_Pi_h[kkk++] = IPJ(2 * e + 1, i, 1);

        Pt_Pi_h[i++] = B * (QFE[p].x) + A * (1. - QFE[p].x);    // X=0 => A  X=1 => B;
      }
    }
    int i6 = 6, i7 = 7;
    for(int p = 0; p < QFK.n; ++p) {
      ipj_Pi_h[kkk++] = IPJ(i6, i, 0);
      ipj_Pi_h[kkk++] = IPJ(i6, i, 1);
      ipj_Pi_h[kkk++] = IPJ(i7, i, 0);
      ipj_Pi_h[kkk++] = IPJ(i7, i, 1);
      Pt_Pi_h[i++] = QFK[p];
    }

    assert(kkk == this->ipj_Pi_h.N());
    assert(i == this->Pt_Pi_h.N());
  }

  void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_<double> &v) const {    // compute the coef of interpolation ...
    const Element &T = K.T;
    int k = 0;
    int arrEdgeOrient[3] = {T.EdgeOrientation(0), T.EdgeOrientation(1), T.EdgeOrientation(2)};
    double s = 1./sqrt(T.measure());
    for(int i = 0; i < 3; i++){
      R2 E(-T.Edge(i).perp());
      R eOrientation = arrEdgeOrient[i];

      for(int p = 0; p < QFE.n; ++p){
        R l0 = QFE[p].x, l1 = 1 - QFE[p].x;
        R p0 = (2 * l0 - l1) * 2;     // poly othogonaux to \lambda_1
        R p1 = (2 * l1 - l0) * 2;     // poly othogonaux to \lambda_0
        R lambda1 = eOrientation * p0 * QFE[p].a;    // [some quadrature function?]
        R lambda0 = eOrientation * p1 * QFE[p].a;    //
        if(eOrientation < 0) {
          Exchange(lambda1, lambda0);    // exch lambda0,lambda1
        }
        v[k++] = lambda0 * E.x * s;
        v[k++] = lambda0 * E.y * s;
        v[k++] = lambda1 * E.x * s;
        v[k++] = lambda1 * E.y * s;
      }
    }

    R2 B[2] = {T.Edge(1), T.Edge(2)};
    B[0] = B[0].perp();
    B[1] = B[1].perp();

    double CK = 0.5;    // dof U= [u1,u2] > |K| int_K ( B_i.U )

    for(int p = 0; p < QFK.n; ++p) {
      double w = QFK[p].a * CK;
      v[k++] = w * B[0].x * s;
      v[k++] = w * B[0].y * s;
      v[k++] = w * B[1].x * s;
      v[k++] = w * B[1].y * s;
    }

    assert(k == this->ipj_Pi_h.N());
  }

  void FB(const What_d, const Element &K, const Rd &PHat, RNMK_ &bfMat) const;


};

// ENDOFCLASS TypedgeOrientationfFE_PkEdge

void TypeOfFE_RT1_2d::FB(const What_d whatd, const Element &K, const Rd &Phat, RNMK_ &bfMat) const{
    R2 X = K(Phat); // [Phat is the quadrature point X mapped to the reference triangle (in loop of addElementMat, X is sent to Phat)]
    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])}; // Triangle K node points
    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y; // [reference triangle barycentric coords]
    R refBaryc[3] = {l0, l1, l2};
    int arrEdgeOrient[] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    /*
    *
    * THE 2 DOF k=0,1  are: on edge e   f -> \int_e f \lambda_{e+k} . n_e
    * THE 2 internal dof are : f -> \int_K f e_i  where e_i is the canonical basis of R^2
    *
    *
    * so the basis function are
    *
    * let call \phi_i the basic fonction of RT0 (without orientation) so the normal is exterior.
    * \phi_i (X) = ( X- Q_i ) / (2 |K|) =  \lambda_{i+1} Curl( \lambda_{i+2}) - \lambda_{i+2} Curl(
    * \lambda_{i+1})
    *
    * edge function j=0,1
    * i1= i+j+1, i2= i+2-j  remark : {i,i1,i2} <=> {i,i+1,i+2}
    * fb_i,j = \phi_i ( \lambda_{i1} - 4/3 \lambda_i) + 1/3 \phi_{i1}\lambda_{i1}
    *
    *
    *
    * we have 2 bubbles functions
    * fb_6 =      8    \phi_0   \lambda_{0} +   16 \phi_1   \lambda_{1}
    * fb_7 =     - 8   \phi_0   \lambda_{0} + 8 \phi_1   \lambda_{1}
    *
    *
    * such that  for i,j =0,1
    * int_K ( Bb_j  f_{6+i))  =   \delta_{ij}
    *
    *
    * so all basic d function are the sum of 3 function
    *
    * sum_{k=0}^2  c_k  phi_{p_k} lambda_{l_k}
    *
    */

    assert(bfMat.N( ) >= ndf);
    assert(bfMat.M( ) == 2);

    bfMat = 0;

    R triMeas2 = 2 * K.mesure();

    R2 phi[3] = {X - Q[0], X - Q[1], X - Q[2]};    // phi * area *2
    R2 X_dx(1,0); R2 X_dy(0,1);
    R2 phi_dx[3] = {X_dx, X_dx, X_dx};
    R2 phi_dy[3] = {X_dy, X_dy, X_dy};

    // Kmap = 1/triMeas2*[Q[2].y-Q[0].y -(Q[2].x-Q[0].x)
    //                  -(Q[1].y-Q[0].y) Q[1].x-Q[0].x]
    R Kmap[2][2]; // Kmap: K -> Khat
    Kmap[0][0] =  1/triMeas2*(Q[2].y-Q[0].y);
    Kmap[0][1] = -1/triMeas2*(Q[2].x-Q[0].x);
    Kmap[1][0] = -1/triMeas2*(Q[1].y-Q[0].y);
    Kmap[1][1] =  1/triMeas2*(Q[1].x-Q[0].x);

    // R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y; // [reference triangle barycentric coords]
    // R refBaryc[3] = {l0, l1, l2};
    R2 Phat_dx( Kmap[0][0], Kmap[1][0] );
    R2 Phat_dy( Kmap[0][1], Kmap[1][1] );
    R refBaryc_dx[3] = {1 - Phat_dx.x - Phat_dx.y, Phat_dx.x, Phat_dx.y};
    R refBaryc_dy[3] = {1 - Phat_dy.x - Phat_dy.y, Phat_dy.x, Phat_dy.y};
    int pI[8][3];    // store p_k
    int lI[8][3];    // store l_k
    R cI[8][3];      // store c_k
    int dof = 0;
    double s = sqrt(K.mesure());

    for (int e = 0; e < 3; ++e) { // [loops through edges]
      //int i = e;
      int ii[2] = {(e + 1) % 3, (e + 2) % 3};
      R eOrientation = arrEdgeOrient[e] / triMeas2;
      if (eOrientation < 0) {
        Exchange(ii[0], ii[1]);
      }

      for (int j = 0; j < 2; ++j, dof++) {
        pI[dof][0] = e;
        lI[dof][0] = ii[j];
        cI[dof][0] = eOrientation;

        pI[dof][1] = e;
        lI[dof][1] = e;
        cI[dof][1] = -eOrientation * 4. / 3.;

        pI[dof][2] = ii[j];
        lI[dof][2] = ii[j];
        cI[dof][2] = eOrientation / 3.;
      }
    } // [now dof = 6, next code piece does bubbles dof=6,7]

    // FB  (x-Q_i) l_i l_j  =
    R s8 = 8 / triMeas2, s01 = s8;
    R cbb[] = {s8, 2 * s01, -s01, s8};    // { [ 8, 16], [ -8, 8] }

    // [the 2 bubbles]
    for (int j = 0; j < 2; ++j, dof++) {    // [j indexes the bubble funs]
      pI[dof][0] = 0;                       // i
      lI[dof][0] = 0;
      cI[dof][0] = cbb[j];

      pI[dof][1] = 1;    // i
      lI[dof][1] = 1;
      cI[dof][1] = cbb[j + 2];

      pI[dof][2] = 2;
      lI[dof][2] = 2;
      cI[dof][2] = 0; // [the third constant is zero for the bubbles]
    }

    assert(dof == 8);

    if (whatd & Fop_id) {
      for (int dof = 0; dof < 8; ++dof) {
        R2 fd(0., 0.);

        for (int k = 0; k < 3; ++k) {
          fd += (cI[dof][k] * refBaryc[lI[dof][k]]) * phi[pI[dof][k]];
        }

        bfMat(dof, 0, op_id) = fd.x * s;
        bfMat(dof, 1, op_id) = fd.y * s;
      }
    }

    if ((whatd & Fop_D1) || (whatd & Fop_D2) ) {
      R2 DL[3] = {K.H(0), K.H(1), K.H(2)}; // [heights?? represents derivative of refBaryc]
      R2 e1(1, 0); // [initialise R2 canonical basis vecs]
      R2 e2(0, 1);

      if (whatd & Fop_dx)  {
        for (int dof = 0; dof < 8; ++dof) {
          R2 fd(0., 0.); // init 0 vec

          for (int k = 0; k < 3; ++k) {
            // std::cout << pI[dof][k] << std::endl;
            fd += cI[dof][k] * (DL[lI[dof][k]].x * phi[pI[dof][k]] + refBaryc[lI[dof][k]] * e1);
          }

          bfMat(dof, 0, op_dx) = fd.x * s;
          bfMat(dof, 1, op_dx) = fd.y * s;
        }
      }

      if (whatd & Fop_dy) {
        for (int dof = 0; dof < 8; ++dof) {
          R2 fd(0., 0.);

          for (int k = 0; k < 3; ++k) {
            fd += cI[dof][k] * (DL[lI[dof][k]].y * phi[pI[dof][k]] + refBaryc[lI[dof][k]] * e2);
          }

          bfMat(dof, 0, op_dy) = fd.x * s;
          bfMat(dof, 1, op_dy) = fd.y * s;
        }
      }

      if (whatd & Fop_D2) {
        if (whatd & Fop_dxx)  {
          for (int dof = 0; dof < 8; ++dof) {
            R2 fd(0., 0.); // init 0 vec

            for (int k = 0; k < 3; ++k) { // [take expression from Fop_dx and diff wrt x]
              // fd += cI[dof][k] * (DL[lI[dof][k]].x * phi_dx[pI[dof][k]] + refBaryc_dx[lI[dof][k]] * e1);
              // fd += cI[dof][k] * (DL[lI[dof][k]].x * phi_dx[pI[dof][k]] + DL[lI[dof][k]].x * e1);
              fd += cI[dof][k] * (DL[lI[dof][k]].x * e1 + DL[lI[dof][k]].x * e1);
            }

            bfMat(dof, 0, op_dxx) = fd.x * s;
            bfMat(dof, 1, op_dxx) = fd.y * s;
          }
        }

        if (whatd & Fop_dxy)  {
          for (int dof = 0; dof < 8; ++dof) {
            R2 fd(0., 0.); // init 0 vec

            for (int k = 0; k < 3; ++k) { // [take expression from Fop_dx and diff wrt y]
              // std::cout << pI[dof][k] << std::endl;
              // fd += cI[dof][k] * (DL[lI[dof][k]].x * phi_dy[pI[dof][k]] + refBaryc_dy[lI[dof][k]] * e1);
              // fd += cI[dof][k] * (DL[lI[dof][k]].x * phi_dy[pI[dof][k]] + DL[lI[dof][k]].y * e1);
              fd += cI[dof][k] * (DL[lI[dof][k]].x * e2 + DL[lI[dof][k]].y * e1);
            }

            bfMat(dof, 0, op_dxy) = fd.x * s;
            bfMat(dof, 1, op_dxy) = fd.y * s;
          }
        }

        if (whatd & Fop_dyy) {
          for (int dof = 0; dof < 8; ++dof) {
            R2 fd(0., 0.);

            for (int k = 0; k < 3; ++k) {
              // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dy[pI[dof][k]] + refBaryc_dy[lI[dof][k]] * e2);
              // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dy[pI[dof][k]] + DL[lI[dof][k]].y * e2);
              fd += cI[dof][k] * (DL[lI[dof][k]].y * e2 + DL[lI[dof][k]].y * e2);
            }

            bfMat(dof, 0, op_dyy) = fd.x * s;
            bfMat(dof, 1, op_dyy) = fd.y * s;
          }
        }

        // if (whatd & Fop_dyx) { // [MIXED PARTIALS DO NOT COMMUTE]
        //   for (int dof = 0; dof < 8; ++dof) {
        //     R2 fd(0., 0.);
        //     // R2 test(0., 0.);
        //
        //     for (int k = 0; k < 3; ++k) {// [take expression from Fop_dy and diff wrt x]
        //       // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dx[pI[dof][k]] + refBaryc_dx[lI[dof][k]] * e2);
        //       // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dx[pI[dof][k]] + DL[lI[dof][k]].x * e2);
        //       fd += cI[dof][k] * (DL[lI[dof][k]].y * e1 + DL[lI[dof][k]].x * e2);
        //     }
        //     // std::cout << (fd-test) << std::endl;
        //     bfMat(dof, 0, op_dyx) = fd.x;
        //     bfMat(dof, 1, op_dyx) = fd.y;
        //   }
        // }
      }
    }
  }

static TypeOfFE_RT1_2d  myRT1_2d;
GTypeOfFE<Mesh2> & RT1_2d(myRT1_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::RT1=myRT1_2d;



class TypeOfFE_RT2_2d : public InitTypeOfRTk_2d, public GTypeOfFE<Mesh2> {
  typedef  Mesh2 Mesh; // Define 2D mesh as Mesh
  typedef  typename Mesh::Element  E; // Define mesh element as E

 public:
  static double Pi_h_coef[];
  bool Ortho;

  TypeOfFE_RT2_2d(bool ortho)
    : InitTypeOfRTk_2d(2), GTypeOfFE<Mesh2>(ndf, 2, Data,
                                    2 * 3 * 3 * QFE.n + QFK.n * 4 * 3,    // nb coef mat interpole
                                    3 * QFE.n + QFK.n,                    // nb P interpolation
                                    0),
      Ortho(ortho) {
    int dofE = this->k + 1;                 // == 3
    int dofKs = (dofE - 1) * (dofE) / 2;    //= = 3 ..(

    ffassert(dofKs == 3);
    ffassert(dofE == 3);

    int kkk = 0, i = 0;

    for (int e = 0; e < 3; ++e) {
      for (int p = 0; p < QFE.n; ++p) {
        R2 A(TriangleHat[Element::nvedge[e][0]]);
        R2 B(TriangleHat[Element::nvedge[e][1]]);

        for (int l = 0; l < dofE; ++l) {
          ipj_Pi_h[kkk++] = IPJ(dofE * e + l, i, 0);
          ipj_Pi_h[kkk++] = IPJ(dofE * e + l, i, 1);
        }

        Pt_Pi_h[i++] = B * (QFE[p].x) + A * (1. - QFE[p].x);    // X=0 => A  X=1 => B;
      }
    }

    for (int p = 0; p < QFK.n; ++p) {
      int i6 = 3 * 3, i7 = i6 + 1;

      for (int l = 0; l < dofKs; ++l) {
        ipj_Pi_h[kkk++] = IPJ(i6, i, 0);
        ipj_Pi_h[kkk++] = IPJ(i6, i, 1);
        ipj_Pi_h[kkk++] = IPJ(i7, i, 0);
        ipj_Pi_h[kkk++] = IPJ(i7, i, 1);
        i6 += 2;
        i7 += 2;
      }

      Pt_Pi_h[i++] = QFK[p];
    }

    assert(kkk == this->ipj_Pi_h.N( ));
    assert(i == this->Pt_Pi_h.N( ));
  }

  void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K,KN_< double > &v) const {
    const Element  &T(K.T);
    int k = 0;
    // magic fom:
    // inv of [[4!,3!,2!2!],[3!,2!2!,2!],[2!2!,3!,4!]]/5! =
    double c1[][3] = {{9, -18, 3} /* 0 */, {-18, 84, -18} /* 1 */, {3, -18, 9} /* 2 */};

    for (int i = 0; i < 3; i++) {
      R2 E(Ortho ? T.Edge(i) : -T.Edge(i).perp( ));
      R s = T.EdgeOrientation(i);

      for (int p = 0; p < QFE.n; ++p) {
        R l1 = QFE[p].x, l2 = 1 - QFE[p].x;
        R l11 = l1 * l1;
        R l22 = l2 * l2;
        R l21 = l2 * l1;
        R p0 = c1[0][0] * l11 + c1[0][1] * l21 + c1[0][2] * l22;    //
        R p1 = c1[1][0] * l11 + c1[1][1] * l21 + c1[1][2] * l22;    //
        R p2 = c1[2][0] * l11 + c1[2][1] * l21 + c1[2][2] * l22;    //
        R sa = s * QFE[p].a;
        R cc2 = sa * p0;    //
        R cc1 = sa * p1;    //
        R cc0 = sa * p2;    //
        if (s < 0) {
          swap(cc0, cc2);
        }

        v[k++] = cc0 * E.x;
        v[k++] = cc0 * E.y;
        v[k++] = cc1 * E.x;
        v[k++] = cc1 * E.y;
        v[k++] = cc2 * E.x;
        v[k++] = cc2 * E.y;
      }
    }

    R2 B[2] = {T.Edge(1), T.Edge(2)};
    if (Ortho) {
      B[0] = -B[0];
      B[1] = -B[1];
    } else {
      B[0] = B[0].perp( );
      B[1] = B[1].perp( );
    }

    double CK = 0.5;    // dof U= [u1,u2] > |K| int_K ( B_i.U )
    R ll[3];            //, lo[3];

    for (int p = 0; p < QFK.n; ++p) {
      double w = -QFK[p].a * CK;
      QFK[p].toBary(ll);
      ll[0] *= w;
      ll[1] *= w;
      ll[2] *= w;

      for (int l = 0; l < 3; ++l) {
        v[k++] = ll[l] * B[0].x;
        v[k++] = ll[l] * B[0].y;
        v[k++] = ll[l] * B[1].x;
        v[k++] = ll[l] * B[1].y;
      }
    }

    assert(k == this->ipj_Pi_h.N( ));
  }

  void FB(const What_d whatd, const Element &K, const Rd &Phat, RNMK_ &val) const;
};

void TypeOfFE_RT2_2d::FB(const What_d whatd, const Element &K, const R2 &Phat, RNMK_ &val) const{
  R2 X = K(Phat);
  R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
  R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y;
  R L[3] = {l0, l1, l2};
  R eo[] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};

  assert(val.N( ) >= ndf);
  assert(val.M( ) == 2);

  val = 0;
  int p[15] = {0, 1, 2,  5,  4,  3,  6, 7,
               8, 9, 10, 11, 12, 13, 14};    // Permutation for orinatation
  R2 Pm[18];                                 // all the momome function ..
  double cf[][6] = {{0, -5.5, 0, -2.5, -0.5, -1.5} /* 0 */,
                    {-1.25, -1.25, 0.25, -1, 0.25, -1} /* 1 */,
                    {-5.5, 0, -0.5, -1.5, 0, -2.5} /* 2 */,
                    {0, -2.5, 0, -5.5, -1.5, -0.5} /* 3 */,
                    {0.25, -1, -1.25, -1.25, -1, 0.25} /* 4 */,
                    {-0.5, -1.5, -5.5, 0, -2.5, 0} /* 5 */,
                    {-2.5, 0, -1.5, -0.5, 0, -5.5} /* 6 */,
                    {-1, 0.25, -1, 0.25, -1.25, -1.25} /* 7 */,
                    {-1.5, -0.5, -2.5, 0, -5.5, 0} /* 8 */,
                    {30, 90, -30, 180, 30, 60} /* 9 */,
                    {90, 30, 30, 60, -30, 180} /* 10 */,
                    {30, -180, -30, -90, -60, -30} /* 11 */,
                    {60, -120, 60, -60, 120, -60} /* 12 */,
                    {-120, 60, 120, -60, 60, -60} /* 13 */,
                    {-180, 30, -60, -30, -30, -90} /* 14 */};
  int Bii[][3] = {{0, 0, 0} /* 0 */,  {0, 1, 1} /* 1 */,  {0, 2, 2} /* 2 */,  {0, 1, 2} /* 3 */,
                  {0, 2, 0} /* 4 */,  {0, 0, 1} /* 5 */,  {1, 0, 0} /* 6 */,  {1, 1, 1} /* 7 */,
                  {1, 2, 2} /* 8 */,  {1, 1, 2} /* 9 */,  {1, 2, 0} /* 10 */, {1, 0, 1} /* 11 */,
                  {2, 0, 0} /* 12 */, {2, 1, 1} /* 13 */, {2, 2, 2} /* 14 */, {2, 1, 2} /* 15 */,
                  {2, 2, 0} /* 16 */, {2, 0, 1} /* 17 */};
  int fe[] = {1, 3, 2, 6, 10, 8, 12, 17, 13};
  int k6[] = {4, 5, 9, 11, 15, 16};
  R CKK = K.mesure() * 2;
  R2 phi[3] = {X - Q[0], X - Q[1], X - Q[2]};    // phi * area *2
  if (Ortho){
    phi[0] = phi[0].perp( );
    phi[1] = phi[1].perp( );
    phi[2] = phi[2].perp( );
  }

  for (int l = 0; l < 18; ++l){
    int i = Bii[l][0];
    int j = Bii[l][1];
    int k = Bii[l][2];
    Pm[l] = phi[i] * (L[j] * L[k] / CKK);
  }

  // static int ddd=0;

  if (eo[0] < 0){
    Exchange(p[0], p[2]);
  }

  if (eo[1] < 0){
    Exchange(p[3], p[5]);
  }

  if (eo[2] < 0){
    Exchange(p[6], p[8]);
  }

  double sg[15] = {eo[0], eo[0], eo[0], eo[1], eo[1], eo[1], eo[2], eo[2],
                   eo[2], 1.,    1.,    1.,    1.,    1.,    1.};

  if (whatd & Fop_id){
    for (int pdf = 0; pdf < 15; ++pdf) {
      int df = p[pdf];
      R2 fd(0., 0.);
      if (df < 9) {
        fd = Pm[fe[df]];    // edge function ..
      }

      for (int k = 0; k < 6; ++k) {
        fd += cf[df][k] * Pm[k6[k]];
      }

      fd *= sg[df];
      val(pdf, 0, op_id) = fd.x;
      val(pdf, 1, op_id) = fd.y;
    }
  }

  if ((whatd & Fop_D1) || (whatd & Fop_D2) ) {
    R2 DL[3] = {K.H(0), K.H(1), K.H(2)};
    R2 Dphi1(1, 0);    // D(phi.x)
    R2 Dphi2(0, 1);    // D(phi.y)
    R2 DxPm[18];
    R2 DyPm[18];

    if (Ortho) {
      Dphi1 = R2(0, -1);
      Dphi2 = R2(1, 0);
    }    // Correction Jan 2019 FH. 	// x,y -> (-y,x)

    for (int l = 0; l < 18; ++l) {
      int i = Bii[l][0];
      int j = Bii[l][1];
      int k = Bii[l][2];
      R Ljk = L[j] * L[k];
      R2 DLjk = L[j] * DL[k] + DL[j] * L[k];
      R2 DF1 = (Dphi1 * Ljk + phi[i].x * DLjk) / CKK;    // BUG ici Ortho ????????
      R2 DF2 = (Dphi2 * Ljk + phi[i].y * DLjk) / CKK;
      DxPm[l] = R2(DF1.x, DF2.x);
      DyPm[l] = R2(DF1.y, DF2.y);
    }

    if (whatd & Fop_dx) {
      for (int pdf = 0; pdf < 15; ++pdf) {
        int df = p[pdf];
        R2 fd(0., 0.);
        if (df < 9) {
          fd = DxPm[fe[df]];    // edge function ..
        }

        for (int k = 0; k < 6; ++k) {
          fd += cf[df][k] * DxPm[k6[k]];
        }

        fd *= sg[df];
        val(pdf, 0, op_dx) = fd.x;
        val(pdf, 1, op_dx) = fd.y;
      }
    }

    if (whatd & Fop_dy) {
      for (int pdf = 0; pdf < 15; ++pdf) {
        int df = p[pdf];
        R2 fd(0., 0.);
        if (df < 9) {
          fd = DyPm[fe[df]];    // edge function ..
        }

        for (int k = 0; k < 6; ++k) {
          fd += cf[df][k] * DyPm[k6[k]];
        }

        fd *= sg[df];
        val(pdf, 0, op_dy) = fd.x;
        val(pdf, 1, op_dy) = fd.y;
      }
    }

    if (whatd & Fop_D2) {
      cout << " to do FH RT2 dxx, dyy dxy " << endl;
      ffassert(0);
    }
  }
}


static TypeOfFE_RT2_2d  myRT2_2d(false);
GTypeOfFE<Mesh2> & RT2_2d(myRT2_2d);
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::RT2=myRT2_2d;
