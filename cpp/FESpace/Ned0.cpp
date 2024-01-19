#include "FESpace.hpp"
#include "transformation.hpp"
// Nedelec type 1 element 3D

class TypeOfFE_Ned0_kind1 : public GTypeOfFE<Mesh3> {
  public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    // typedef GFElement<Mesh3> FElement;

    // static int dfon[];
    static const int d = Mesh::Rd::d;
    const QuadratureFormular1d &QFE; //   static const GQuadratureFormular<R1> QFe;
    TypeOfFE_Ned0_kind1();
    // int edgeface[4][3] ;
    static int Data[];

    void FB(const What_d whatd, const Element &K, const R3 &PHat, RNMK_ &bfMat) const;
    void get_Coef_Pi_h(const GbaseFElement<Mesh> &K, KN_<double> &v) const; // ,int ocoef,int odf,int *nump  ) const;
};

// int TypeOfFE_Ned0_kind1::dfon[]={0,1,0,0};
int TypeOfFE_Ned0_kind1::Data[]{
    // geometry=vertex,edge,face,volume; 0,1,2 vertices, 3,4,5 edges, 6 face
    4, 5, 6, 7, 8, 9, // on what hard coded geometry the dofs are located
    0, 0, 0, 0, 0, 0, // the (number of the dof - 1) on corresponding geometry/node (here each edge has 1 dof)
    0, 1, 2, 3, 4, 5, // the geometry(vertex) of the basis fcn corresponding to the dof
    0, 1, 2, 3, 4, 5, // the dof of the sub FE (default 0,1,2) [???]
    0, 1, 0, 0,       // which geometry=vertex,edge,face,volume contains dofs
    0,                // for each component $j=0,N-1$ it give the sub FE associated [???]
    0,                // begin_dfcomp
    6                 // end_dfcomp (=total number of dof in one element)
};
// const  GQuadratureFormular<R1> TypeOfFE_Ned0_kind1::QFe(-1+2*2,2,GaussLegendre(2),true);

TypeOfFE_Ned0_kind1::TypeOfFE_Ned0_kind1()
    : GTypeOfFE<Mesh3>(6, 3, Data,
                       18 * 2, // Element::ne*3*QFe.n,
                       6 * 2), // Element::ne*QFe.n)
      QFE(QF_GaussLegendre2) {
    // assert(QFe.n);

    GTypeOfFE<Mesh>::basisFctType    = BasisFctType::Ned0;
    GTypeOfFE<Mesh>::polynomialOrder = 1;

    //  integration on edge use QFe
    R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.)};
    int i   = 0;
    int kkk = 0;
    for (int e = 0; e < 6; ++e) { // loop over edges
        for (int q = 0; q < QFE.n; ++q, ++i) {
            double x   = QFE[q].x;
            Pt_Pi_h[i] = Pt[Element::nvedge[e][0]] * x + Pt[Element::nvedge[e][1]] * (1 - x);

            ipj_Pi_h[kkk++] = IPJ(e, i, 0);
            ipj_Pi_h[kkk++] = IPJ(e, i, 1);
            ipj_Pi_h[kkk++] = IPJ(e, i, 2);
        }
    }
    assert(kkk == 36);
    assert(i == 12);
    // int i=0,p=0;
    // for (int e=0;e<6;e++) {
    //   for(int q=0;q<QFE.n;++q,++p) {
    //     for (int c=0;c<3;c++,i++) {
    //       ipj_Pi_h[i] = IPJ(e, i, c);

    //       this->pInterpolation[i]=p;
    //       this->cInterpolation[i]=c;
    //       this->dofInterpolation[i]=e;
    //       this->coefInterpolation[i]=0.;

    //       ipj_Pi_h[i] = IPJ(p,e,c);
    //     }
    //   }
    // }
    // cout <<  " ++ TypeOfFE_Ned0_kind1():"<< this->PtInterpolation << endl;
}

void TypeOfFE_Ned0_kind1::get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                                        KN_<double> &v) const // ,int ocoef,int odf,int *nump ) const
{
    // compute de coef d'interpolation
    const Element &T = K.T;
    int i = 0, p = 0; // int i=ocoef,p=0;
    for (int e = 0; e < 6; ++e) {

        //! Added
        // int idx_face = T.faceOfEdge[e][0];
        // if (T.EdgeOrientation(e) < 0) {     
        //     idx_face = T.faceOfEdge[e][1];
        // }
        // R cc = T.N_notNormalized(idx_face).norm(); // area of face 0
        
        R3 E = T.EdgeOrientation(e) * T.Edge(e); //  exterior and  ||N|| = 2* area f
        //R cc = T.Edge(e).norm();
        R cc = 1.;

        for (int q = 0; q < QFE.n; ++q, ++p) {
            for (int c = 0; c < 3; c++, i++) {
                //v[i] = E[c] * QFE[q].a;   //! Original

                int idx_face = T.faceOfEdge[i][0];

                //v[i] = E[c] * QFE[q].a / (T.N_notNormalized(0).norm() / 2);      //! Mine
                v[i] = E[c] * QFE[q].a / cc;      //! Mine
                
                
            }
        }
        // ffassert(i==M.ncoef && M.np == p );
    }
}

void TypeOfFE_Ned0_kind1::FB(const What_d whatd, const Element &K, const R3 &PHat, RNMK_ &bfMat) const {
    assert(bfMat.N() >= 6);
    assert(bfMat.M() == 3);
    R l[] = {1. - PHat.sum(), PHat.x, PHat.y, PHat.z};
    R3 D[4];
    K.Gradlambda(D);    // this seems to divide by 6*volume(K)  (6 = 2*dim)

    // wi = signe * (x - qi) / (volume*d)
    bfMat = 0;
    //  i,j : l1 grad lj - lj grad lj
    // int_i^j  grad lj . t_ij = 1

    int se[] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2),
                K.EdgeOrientation(3), K.EdgeOrientation(4), K.EdgeOrientation(5)};

    //R cc = K.N_notNormalized(0).norm() / 2; // area of face 0
    R cc = 1.;

    if (whatd & Fop_D0) {
        R3 X = K(PHat);
        // int k=0;
        for (int i = 0; i < 6; ++i) {
            int i0 = Element::nvedge[i][0], i1 = Element::nvedge[i][1];
            if (se[i] < 0) 
                std::swap(i0, i1);   //! Original
            
            // int idx_face = K.faceOfEdge[i][0];
            // if (se[i] < 0) {    //! Mine
            //     std::swap(i0, i1);   
            //     idx_face = K.faceOfEdge[i][1];
            // }
            
            //R3 wi              = l[i0] * D[i1] - l[i1] * D[i0];   //! Original
            
            //R cc = K.N_notNormalized(idx_face).norm(); // area of face 0 (should not divide by 2 since we divide by 2 in K.gradLambda)
            // std::cout << "cc in FB: " << cc << "\n";
            // std::cout << "face 0: " << K.N_notNormalized(0).norm() << "\n";
            // R3 wi              = (l[i0] * D[i1] - l[i1] * D[i0]) * cc; //! Mine
            R3 wi              = (l[i0] * D[i1] - l[i1] * D[i0]) * K.Edge(i).norm(); //! Mine
            bfMat(i, 0, op_id) = wi.x;
            bfMat(i, 1, op_id) = wi.y;
            bfMat(i, 2, op_id) = wi.z;
            // cout  << "Edge0 3d  "<<i << " "<< X << " " <<wi << " fo: " <<se[i] <<endl;
        }
    }

    if (whatd & Fop_D1) {
        for (int i = 0; i < 6; ++i) {
            int i0 = Element::nvedge[i][0], i1 = Element::nvedge[i][1];
            
            //! Original
            // if (se[i] < 0)
            //     std::swap(i0, i1);

            //! Mine
            int idx_face = K.faceOfEdge[i][0];
            if (se[i] < 0) {    
                std::swap(i0, i1);   
                idx_face = K.faceOfEdge[i][1];
            }
            //R cc = K.N_notNormalized(idx_face).norm(); // area of face 0
            R cc = K.Edge(i).norm();

            if (whatd & Fop_dx) {
                //R3 wi              = D[i0].x * D[i1] - D[i1].x * D[i0];   //! Original
                R3 wi              = (D[i0].x * D[i1] - D[i1].x * D[i0]) * cc; // ! Mine
                bfMat(i, 0, op_dx) = wi.x;
                bfMat(i, 1, op_dx) = wi.y;
                bfMat(i, 2, op_dx) = wi.z;
            }
            if (whatd & Fop_dy) {
                //R3 wi              = D[i0].y * D[i1] - D[i1].y * D[i0];  //! Original
                R3 wi              = (D[i0].y * D[i1] - D[i1].y * D[i0]) * cc; //! Mine
                bfMat(i, 0, op_dy) = wi.x;
                bfMat(i, 1, op_dy) = wi.y;
                bfMat(i, 2, op_dy) = wi.z;
            }
            if (whatd & Fop_dz) {
                //R3 wi              = D[i0].z * D[i1] - D[i1].z * D[i0];  //! Original
                R3 wi              = (D[i0].z * D[i1] - D[i1].z * D[i0]) * cc; //! Mine
                bfMat(i, 0, op_dz) = wi.x;
                bfMat(i, 1, op_dz) = wi.y;
                bfMat(i, 2, op_dz) = wi.z;
            }
        }
    }
}

static TypeOfFE_Ned0_kind1 Ned0_kind1;
GTypeOfFE<Mesh3> &Ned0kind1(Ned0_kind1);
template <> GTypeOfFE<Mesh3> &DataFE<Mesh3>::Ned0 = Ned0kind1;