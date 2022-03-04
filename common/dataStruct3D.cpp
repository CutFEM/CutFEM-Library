#include "dataStruct3D.hpp"

const R3 R3::KHat[4]={R3(0.,0.,0.),R3(1.,0.,0.),R3(0.,1.,0.),R3(0.,0.,1.)};

  //  Attention  nvfaceTet  donnne les faces  les 4 faces de tet telle que la
  // tel que  le tet forme des trois sommet  + l'autre sommet soit positif.
  //  donc  le  produit vectoriel des 2 aretes  (0,1) (0,2)  donne une  normale interieur.
  //  Ok pour les gradients des $\lambda_i$
static const int  nvfaceTet[4][3] ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}   ;
static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };

static const int  nvfaceTria[1][3] = { {0,1,2} };
static const int  nvedgeTria[3][2] = { {1,2},{2,0},{0,1}};

static const int  nvfaceSeg[1][3]  = {{-1,-1,1}};
static const int  nvedgeSeg[1][2]  = { {0,1} };

static const int  nvfaceQuad3[1][4]  = { {0,1,2,3} };
static const int  nvedgeQuad3[4][2] = {{0,1},{1,2},{2,3},{3,0}};

static const int  nvfaceHexa[6][4] ={{0,1,2,3}, {0,1,5,4}, {1,2,6,5},
                                     {2,3,7,6}, {3,0,4,7}, {4,5,6,7}}   ;
static const int  nvedgeHexa[12][2] = { {0,1},{1,2},{2,3},{3,0},
                                        {0,4},{1,5},{2,6},{3,7},
                                        {4,5},{5,6},{6,7},{7,4}};
static const int commonVertOfEdgeHexa[12][12] = {
   {-1,1,-1,0,0,1,-1,-1,-1,-1,-1,-1},
   {1,-1,2,-1,-1,1,2,-1,-1,-1,-1,-1},
   {-1,2,-1,3,-1,-1,2,3,-1,-1,-1,-1},
   {0,-1,3,-1,0,-1,-1,3,-1,-1,-1,-1},
   {0,-1,-1,0,-1,-1,-1,-1,4,-1,-1,4},
   {1,1,-1,-1,-1,-1,-1,-1,5,5,-1,-1},
   {-1,2,2,-1,-1,-1,-1,-1,-1,6,6,-1},
   {-1,-1,3,3,-1,-1,-1,-1,-1,-1,7,7},
   {-1,-1,-1,-1,4,5,-1,-1,-1,5,-1,4},
   {-1,-1,-1,-1,-1,5,6,-1,5,-1,6,-1},
   {-1,-1,-1,-1,-1,-1,6,7,-1,6,-1,7},
   {-1,-1,-1,-1,4,-1,-1,7,4,-1,7,-1}
   };

static const int faceOfEdgeHexa[12][2] = {{0,1},{0,2},{0,3},{0,4},
                                          {1,4}, {1,2},{2,3},{3,4},
                                          {1,5}, {2,5},{3,5},{4,5}};
static const int edgeOfFaceHexa[6][4]  = {{0,1,2,3}, {0,5,8,4}, {1,6,9,5},
                                          {2,7,10,6}, {3,4,11,7}, 8,9,10,11  };
static const int connectivityNodeHexa[8][3] = {
  {1,3,4},{0,2,5},{1,3,6},{0,2,7},{0,5,7},{1,4,6},{2,5,7},{3,4,6}
} ;

static const int commonVertOfEdgeTria[3][3] = {{-1,2,1} , {2,-1,0} , {1,0,-1}};
static const int commonVertOfEdgeQuad3[4][4] = { {-1,1,-1,0} , {1,-1,2,-1} , {-1,2,-1,3}, {0,-1,3,-1}};
static const int commonVertOfEdgeTet[6][6] = {
  {-1,0,0,1,1,-1} , {0,-1,0,2,-1,2} , {0,0,-1,-1,3,3},
  {1,2,-1,-1,1,2} , {1,-1,3,1,-1,3} , {-1,2,3,2,3,-1}};



static const int nvrefTria[4][3] = { {0,5,4}, {1,3,5}, {2,4,3}, {3,4,5}};

static const int nvrefTet[8][4] = { {0,4,5,6}, {1,7,4,8}, {2,5,7,9}, {3,6,9,8},
				 {4,7,5,8}, {5,8,6,4}, {5,8,9,6}, {5,7,9,8}};

static const int onWhatIsEdge3[3][7] = {
  {0,1,3, 2,0,0, 0}, // edge 0
  {3,0,1, 0,2,0, 0},
  {1,3,0, 0,0,2, 0}
};
static const int onWhatIsFace3[4][15] = {
  {0,1,1,1, 0,0,0,2,2,2, 3,0,0,0, 0 },
  {1,0,1,1, 0,2,2,0,0,2, 0,3,0,0, 0 },
  {1,1,0,1, 2,0,2,0,2,0, 0,0,3,0, 0 },
  {1,1,1,0, 2,2,0,2,0,0, 0,0,0,3, 0 }
};
static const int onWhatIsEdgeQuad3[4][9] = {
  {1,1,0,0, 2,0,0,0, 0},
  {0,1,1,0, 0,2,0,0, 0},
  {0,0,1,1, 0,0,2,0, 0},
  {1,0,0,1, 0,0,0,2, 0}
};



template<> const int (* const GenericElement<DataTriangle3>::nvface)[3] = nvfaceTria ;
template<> const int (* const GenericElement<DataTriangle3>::nvedge)[2] = nvedgeTria ;
template<> const int (* const GenericElement<DataTriangle3>::nvadj)[2] = nvedgeTria ;
template<> const int  GenericElement<DataTriangle3>::nitemdim[4] = {3,3,1,0 }  ;
template<> const int (* const GenericElement<DataTriangle3>::onWhatBorder)[7] = onWhatIsEdge3 ;
template<> const int (* const GenericElement<DataTriangle3>::commonVertOfEdges)[3] = commonVertOfEdgeTria ;


template<> const int (* const GenericElement<DataQuad3>::nvface)[4] = nvfaceQuad3 ;
template<> const int (* const GenericElement<DataQuad3>::nvedge)[2] = nvedgeQuad3 ;
template<> const int (* const GenericElement<DataQuad3>::nvadj)[2]  = nvedgeQuad3 ;
template<> const int  GenericElement<DataQuad3>::nitemdim[4] = {4,4,1,0 }  ;
template<> const int (* const GenericElement<DataQuad3>::commonVertOfEdges)[4] = commonVertOfEdgeQuad3 ;
template<> const int (* const GenericElement<DataQuad3>::onWhatBorder)[9] = onWhatIsEdgeQuad3 ;


template<> const int (* const GenericElement<DataTet>::nvface)[3] = nvfaceTet ;
template<> const int (* const GenericElement<DataTet>::nvedge)[2] = nvedgeTet ;
template<> const int (* const GenericElement<DataTet>::nvadj)[3] = nvfaceTet ;
template<> const int  GenericElement<DataTet>::nitemdim[4] = {4,6,4,1 }  ;
template<> const int (* const GenericElement<DataTet>::onWhatBorder)[15] = onWhatIsFace3 ;
template<> const int (* const GenericElement<DataTet>::commonVertOfEdges)[6] = commonVertOfEdgeTet ;
const int  Tet::oppEdgeOfEdge[6] = {5,4,3,2,1,0}  ;


template<> const int (* const GenericElement<DataHexa>::nvedge)[2] = nvedgeHexa ;
template<> const int (* const GenericElement<DataHexa>::nvface)[4] = nvfaceHexa;
template<> const int (* const GenericElement<DataHexa>::nvadj)[4] = nvfaceHexa;
// template<> const int (* const GenericElement<DataHexa>::onWhatBorder)[27] = ;
template<> const int (* const GenericElement<DataHexa>::commonVertOfEdges)[12] = commonVertOfEdgeHexa;
template<> const int (* const GenericElement<DataHexa>::faceOfEdge)[2] = faceOfEdgeHexa;
template<> const int (* const GenericElement<DataHexa>::edgeOfFace)[4] = edgeOfFaceHexa;
template<> const int  GenericElement<DataHexa>::nitemdim[4] = {8,12,6,1 }  ;
const int  Hexa::oppEdgeOfEdge[12] = {10,11,8,9,6,7,4,5,2,3,0,1};
const int (*const Hexa::nodeConnectivity)[3] = connectivityNodeHexa;



R3 Tet::H(int i) const
{ ASSERTION(i>=0 && i <4);
  R3 AB(at(this->nvface[i][0]),at(this->nvface[i][1]));
  R3 AC(at(this->nvface[i][0]),at(this->nvface[i][2]));
  return AB^AC/(6.*this->mesure());} // heigth

R3 Tet::n(int i) const
{ ASSERTION(i>=0 && i <4);
  R3 AB(at(this->nvface[i][0]),at(this->nvface[i][1]));
  R3 AC(at(this->nvface[i][0]),at(this->nvface[i][2]));
  R3 N=AB^AC;
  return N/N.norme();} //  exterior normal


R3 Tet::toKref(const R2& P, int i) const {
    return (1 - P.x - P.y) * R3::KHat[nvface[i][0]] + P.x * R3::KHat[nvface[i][1]]
      + P.y * R3::KHat[nvface[i][2]];
  }

R Tet::mesureBord(int i) const {
  R3 AB(at(nvface[i][0]),at(nvface[i][1]));
  R3 AC(at(nvface[i][0]),at(nvface[i][2]));
  return (AB^AC).norm()*0.5;
}
