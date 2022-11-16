#ifndef _GTYPE_OF_FE_HPP
#define _GTYPE_OF_FE_HPP

#include <iostream>
#include <cassert>

#include "../util/util.hpp"
#include "../common/RNM.hpp"
#include "../common/Mesh1dn.hpp"
#include "operationTypeFE.hpp"
// #include "transformation.hpp"

enum class BasisFctType {
   P0,
   P1,
   P2,
   P3,
   P4,
   P1dc,
   P2dc,
   P3dc,
   P0poly,
   P1poly,
   P2poly,
   P2BR,
   RT0,
   RT1,
   RT2,
   BDM1,
   UNDEFINED
};

/*
 *   Basis class for the FE
 *   Contains divers information
 */
class dataTypeOfFE {

 private:
   const int *data; // array to save all the data needed for the FE
   const int *dataalloc;

 public:
   const int ndfonVertex; // number of dof on points
   const int ndfonEdge;   // number of dof on egdes
   const int ndfonFace;   // number of dof on faces
   const int ndfonVolume; // number of dof in volume
   const int nbDoF;       // dof per element
   const int nbNode;      // node per element = dof so far
   int N;                 // dim espace arrivée N = 1
   int nbOfFE;
   bool constDfPerNode = true;

   int const *const nbNodeOnWhat;
   int const *const DFOnWhat;
   int const *const DFOfNode;
   int const *const NodeOfDF;

   const int *ndfOn() const { return &ndfonVertex; }

   dataTypeOfFE(const std::vector<int> &nitemdim,
                const KN<dataTypeOfFE const *> &);
   dataTypeOfFE(const std::vector<int> &nitem, const int *Data, int nbdf,
                int NN);

   virtual ~dataTypeOfFE() {
      if (dataalloc)
         delete[] dataalloc;
   }
};

template <class RdHat> class InterpolationMatrix;
template <class Mesh> class GbaseFElement;
template <class Mesh> class GFElement;

struct IPJ {
   int i, p, j; // i is DoF, p is Point, j is componante
   IPJ(int ii = 0, int pp = 0, int jj = 0) : i(ii), p(pp), j(jj) {}
};

/*
 *   Main class for the FE structure
 *
 */
template <class Mesh> class GTypeOfFE : public dataTypeOfFE {
 public:
   typedef typename Mesh::Element Element;
   typedef Element E;
   typedef typename Element::RdHat RdHat;
   typedef typename Element::Rd Rd;
   typedef GFElement<Mesh> FElement;

   BasisFctType basisFctType = BasisFctType::UNDEFINED;
   int polynomialOrder;

   int NbPtforInterpolation;   // Nb of interpolation points per elemen
   int NbcoefforInterpolation; // Nb of interpolation points per element

   KN<IPJ> ipj_Pi_h;
   KN<RdHat> Pt_Pi_h;
   double *coef_Pi_h;
   // KN<int> begin_coef_Pi_h, end_coef_Pi_h;

   KN<GTypeOfFE<Mesh> *> Sub_ToFE;
   KN<int> begin_dfcomp, end_dfcomp;

   virtual void get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                              KN_<double> &v) const {
      assert(coef_Pi_h);
      // std::cout << "get coef in GTypeOfFE" << std::endl;
      v = KN_<double>(coef_Pi_h, ipj_Pi_h.N());
   }

   /*
    * Simple constructor for simple finite element
    * (not for union of FE)
    */
   GTypeOfFE(const int nbdf, const int NN, const int *data, int kPi, int npPi,
             double *coef_Pi_h_a = 0)
       : dataTypeOfFE(Element::itemTopology(), data, nbdf, NN),

         NbPtforInterpolation(npPi), NbcoefforInterpolation(kPi), ipj_Pi_h(kPi),
         Pt_Pi_h(npPi), coef_Pi_h(coef_Pi_h_a),
         // begin_coef_Pi_h, end_coef_Pi_h;
         Sub_ToFE(nbOfFE), begin_dfcomp(N, 0), end_dfcomp(N, nbdf)
   // begin_dfcomp(data+4*nbdf+4+N),
   // end_dfcomp(data+4*nbdf+4+2*N)
   {
      Sub_ToFE = this;
   }

   /*
    *  Constructor for GTypeOfFE_Sum
    *
    */
   GTypeOfFE(const KN<GTypeOfFE<Mesh> const *> &t)
       : dataTypeOfFE(Element::itemTopology(), t),
         NbPtforInterpolation(this->nbNode),
         NbcoefforInterpolation(this->nbNode), Sub_ToFE(nbOfFE),
         begin_dfcomp(N, 0), end_dfcomp(N, this->nbDoF)

   {

      Sub_ToFE = this;
   }

   //   /*
   //  *  Constructor for GTypeOfFE_Time
   //  *
   //  */
   // GTypeOfFE(const KN<GTypeOfFE<Mesh> const *> &t, const GTypeOfFE<Mesh1>*
   // tt)
   //   :
   //   dataTypeOfFE(Element::nitemdim,t, tt),
   //   NbPtforInterpolation(this->nbNode),
   //   NbcoefforInterpolation(this->nbNode),
   //   PtInterpolation(0),
   //   coefInterpolation(0,0),
   //   Sub_ToFE(nbOfFE),
   //   begin_dfcomp(N,0),
   //   end_dfcomp(N,this->nbDoF)

   // {
   //   Sub_ToFE = this;
   // }

   // virtual void init(InterpolationMatrix<RdHat> & M) const;
   virtual void FB(const What_d whatd, const Element &K, const Rd &P,
                   KNMK_<R> &val) const = 0;
   virtual void FB(const What_d whatd, const Element &K, const Rd &P,
                   KNMK_<R> &val, const KNM_<R> &J) const {
      assert(0);
   };

   ~GTypeOfFE() {}

 private:
   GTypeOfFE(const GTypeOfFE &);
   void operator=(const GTypeOfFE &);
};

// template<class Mesh>
// void GTypeOfFE<Mesh>::init(InterpolationMatrix<RdHat> & M) const
//   {
//     assert(M.np==NbPtforInterpolation);
//     assert(M.ncoef==NbcoefforInterpolation);
//
//     M.P=PtInterpolation;
//     M.coef=coefInterpolation;
//   }

/*
 *   Structure that will contain the different FE
 */
template <class mesh> struct DataFE {
   static GTypeOfFE<mesh> &P0;
   static GTypeOfFE<mesh> &P1;
   static GTypeOfFE<mesh> &P2;
   static GTypeOfFE<mesh> &P3;
   static GTypeOfFE<mesh> &P4;

   static GTypeOfFE<mesh> &P0Poly;
   static GTypeOfFE<mesh> &P1Poly;
   static GTypeOfFE<mesh> &P2Poly;
   static GTypeOfFE<mesh> &P3Poly;

   static GTypeOfFE<mesh> &RT0;
   static GTypeOfFE<mesh> &RT0m;
   static GTypeOfFE<mesh> &RT1;
   static GTypeOfFE<mesh> &RT2;
   static GTypeOfFE<mesh> &BDM1;
   static GTypeOfFE<mesh> &P2BR;

   static GTypeOfFE<mesh> &P1dc;
   static GTypeOfFE<mesh> &P2dc;
   static GTypeOfFE<mesh> &P3dc;

   static GTypeOfFE<mesh> &P0sc;
   static GTypeOfFE<mesh> &P1dcsc;

   static GTypeOfFE<mesh> &P1dcTaylor;
};

#endif
