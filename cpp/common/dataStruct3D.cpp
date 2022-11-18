#include "dataStruct3D.hpp"

const std::vector<R3> R3::KHat = {R3(0., 0., 0.), R3(1., 0., 0.),
                                  R3(0., 1., 0.), R3(0., 0., 1.)};

/// @brief Static mumber for Point3
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint3>::nvedge = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint3>::nvface = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint3>::nvhyperFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint3>::edgeOfFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint3>::faceOfEdge = {
    {}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataPoint3>::commonVertOfEdges = {{}};

/// @brief Static mumber for Seg3
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg3>::nvedge = {{0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg3>::nvface = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg3>::nvhyperFace = {
    {0}, {1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg3>::edgeOfFace = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg3>::faceOfEdge = {{}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataSeg3>::commonVertOfEdges = {{}};

/// @brief Static mumber for Triangle3
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle3>::nvedge = {
    {1, 2}, {2, 0}, {0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle3>::nvface = {
    {0, 1, 2}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle3>::nvhyperFace =
    {{1, 2}, {2, 0}, {0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle3>::edgeOfFace =
    {{0, 1, 2}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle3>::faceOfEdge{
    {0}, {0}, {0}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataTriangle3>::commonVertOfEdges = {
        {-1, 2, 1}, {2, -1, 0}, {1, 0, -1}};

/// @brief Static mumber for Quad3
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad3>::nvedge = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad3>::nvface = {
    {0, 1, 2, 3}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad3>::nvhyperFace = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad3>::edgeOfFace = {
    {0, 1, 2, 3}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad3>::faceOfEdge{
    {0}, {0}, {0}, {0}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataQuad3>::commonVertOfEdges = {
        {-1, 1, -1, 0}, {1, -1, 2, -1}, {-1, 2, -1, 3}, {0, -1, 3, -1}};

/// @brief Static mumber for Tet
template <>
const std::vector<std::vector<int>> GenericElement<DataTet>::nvedge = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTet>::nvface = {
    {3, 2, 1}, {0, 2, 3}, {3, 1, 0}, {0, 1, 2}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTet>::nvhyperFace = {
    {3, 2, 1}, {0, 2, 3}, {3, 1, 0}, {0, 1, 2}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTet>::edgeOfFace = {
    {3, 4, 5}, {1, 2, 5}, {0, 2, 4}, {0, 1, 3}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTet>::faceOfEdge = {
    {2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTet>::commonVertOfEdges =
    {{-1, 0, 0, 1, 1, -1}, {0, -1, 0, 2, -1, 2}, {0, 0, -1, -1, 3, 3},
     {1, 2, -1, -1, 1, 2}, {1, -1, 3, 1, -1, 3}, {-1, 2, 3, 2, 3, -1}};

/// @brief Static mumber for Hexa
template <>
const std::vector<std::vector<int>> GenericElement<DataHexa>::nvedge = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
    {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
template <>
const std::vector<std::vector<int>> GenericElement<DataHexa>::nvface = {
    {0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
    {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
template <>
const std::vector<std::vector<int>> GenericElement<DataHexa>::nvhyperFace = {
    {0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
    {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
template <>
const std::vector<std::vector<int>> GenericElement<DataHexa>::edgeOfFace = {
    {0, 1, 2, 3},  {0, 5, 8, 4},  {1, 6, 9, 5},
    {2, 7, 10, 6}, {3, 4, 11, 7}, {8, 9, 10, 11}};
template <>
const std::vector<std::vector<int>> GenericElement<DataHexa>::faceOfEdge = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 4}, {1, 2},
    {2, 3}, {3, 4}, {1, 5}, {2, 5}, {3, 5}, {4, 5}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataHexa>::commonVertOfEdges = {
        {-1, 1, -1, 0, 0, 1, -1, -1, -1, -1, -1, -1},
        {1, -1, 2, -1, -1, 1, 2, -1, -1, -1, -1, -1},
        {-1, 2, -1, 3, -1, -1, 2, 3, -1, -1, -1, -1},
        {0, -1, 3, -1, 0, -1, -1, 3, -1, -1, -1, -1},
        {0, -1, -1, 0, -1, -1, -1, -1, 4, -1, -1, 4},
        {1, 1, -1, -1, -1, -1, -1, -1, 5, 5, -1, -1},
        {-1, 2, 2, -1, -1, -1, -1, -1, -1, 6, 6, -1},
        {-1, -1, 3, 3, -1, -1, -1, -1, -1, -1, 7, 7},
        {-1, -1, -1, -1, 4, 5, -1, -1, -1, 5, -1, 4},
        {-1, -1, -1, -1, -1, 5, 6, -1, 5, -1, 6, -1},
        {-1, -1, -1, -1, -1, -1, 6, 7, -1, 6, -1, 7},
        {-1, -1, -1, -1, 4, -1, -1, 7, 4, -1, 7, -1}};

const std::vector<int> Tet::oppEdgeOfEdge  = {5, 4, 3, 2, 1, 0};
const std::vector<int> Hexa::oppEdgeOfEdge = {10, 11, 8, 9, 6, 7,
                                              4,  5,  2, 3, 0, 1};
const std::vector<std::vector<int>> Hexa::nodeConnectivity = {
    {1, 3, 4}, {0, 2, 5}, {1, 3, 6}, {0, 2, 7},
    {0, 5, 7}, {1, 4, 6}, {2, 5, 7}, {3, 4, 6}};

R3 Tet::H(int i) const {
   assert(i >= 0 && i < 4);
   R3 AB(at(this->nvface[i][0]), at(this->nvface[i][1]));
   R3 AC(at(this->nvface[i][0]), at(this->nvface[i][2]));
   return AB ^ AC / (6. * this->mesure());
} // heigth

R3 Tet::n(int i) const {
   assert(i >= 0 && i < 4);
   R3 AB(at(this->nvface[i][0]), at(this->nvface[i][1]));
   R3 AC(at(this->nvface[i][0]), at(this->nvface[i][2]));
   R3 N = AB ^ AC;
   return N / N.norme();
} //  exterior normal

R3 Tet::toKref(const R2 &P, int i) const {
   return (1 - P.x - P.y) * R3::KHat[nvface[i][0]] +
          P.x * R3::KHat[nvface[i][1]] + P.y * R3::KHat[nvface[i][2]];
}
R3 Tet::mapToReferenceElement(const R2 &P, int i) const {
   return (1 - P.x - P.y) * R3::KHat[nvface[i][0]] +
          P.x * R3::KHat[nvface[i][1]] + P.y * R3::KHat[nvface[i][2]];
}

R Tet::mesureBord(int i) const {
   R3 AB(at(nvface[i][0]), at(nvface[i][1]));
   R3 AC(at(nvface[i][0]), at(nvface[i][2]));
   return (AB ^ AC).norm() * 0.5;
}

int Tet::faceOrient(int i) const { // def the permutatution of orient the face
   int fo             = 1;
   const Vertex *f[3] = {&at(nvface[i][0]), &at(nvface[i][1]),
                         &at(nvface[i][2])};
   if (f[0] > f[1])
      fo = -fo, std::swap(f[0], f[1]);
   if (f[1] > f[2]) {
      fo = -fo, std::swap(f[1], f[2]);
      if (f[0] > f[1])
         fo = -fo, std::swap(f[0], f[1]);
   }
   return fo;
}
