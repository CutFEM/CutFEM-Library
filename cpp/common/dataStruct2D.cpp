#include "dataStruct2D.hpp"

const std::vector<R2> R2::KHat = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};

static const int nvfaceTet[4][3] = {{2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
static const int nvedgeTet[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {0, 1}, {1, 2}, {2, 3}};

static const int nvfaceTria[1][3] = {{0, 1, 2}};
static const int nvedgeTria[3][2] = {{1, 2}, {2, 0}, {0, 1}};

static const int nvfaceQuad[1][4] = {{0, 1, 2, 3}};
static const int nvedgeQuad[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

static const int nvfaceSeg[1][3] = {{-1, -1, 1}};
static const int nvedgeSeg[1][2] = {{0, 1}};

static const int nvadjSeg[2][1] = {{0}, {1}};

static const int nvref[4][3]    = {{0, 5, 4}, {1, 3, 5}, {2, 4, 3}, {3, 4, 5}};
static const int nvrefSeg[2][2] = {{0, 1}, {1, 2}};

static const int commonVertOfEdge[3][3] = {{-1, 2, 1}, {2, -1, 0}, {1, 0, -1}};
static const int onWhatIsEdge[3][7]     = {{0, 1, 3, 2, 0, 0, 0}, // edge 0
                                           {3, 0, 1, 0, 2, 0, 0},
                                           {1, 3, 0, 0, 0, 2, 0}};
static const int commonVertOfEdgeQuad[4][4] = {
    {-1, 1, -1, 0}, {1, -1, 2, -1}, {-1, 2, -1, 3}, {0, -1, 3, -1}};
static const int onWhatIsEdgeQuad[4][9] = {
    {1, 1, 0, 0, 2, 0, 0, 0, 0}, // edge 0
    {0, 1, 1, 0, 0, 2, 0, 0, 0},
    {0, 0, 1, 1, 0, 0, 2, 0, 0},
    {1, 0, 0, 1, 0, 0, 0, 2, 0}};

template <>
const int (*const GenericElement<DataTriangle2>::nvface)[3] = nvfaceTria;
template <>
const int (*const GenericElement<DataTriangle2>::nvedge)[2] = nvedgeTria;
template <>
const int (*const GenericElement<DataTriangle2>::nvhyperFace)[2] = nvedgeTria;
template <>
const int (*const GenericElement<DataTriangle2>::nvadj)[2]       = nvedgeTria;
template <> const int GenericElement<DataTriangle2>::nitemdim[4] = {3, 3, 1, 0};
template <>
const int (*const GenericElement<DataTriangle2>::commonVertOfEdges)[3] =
    commonVertOfEdge;
// template<> const int (* const GenericElement<DataTriangle2>::refElement)[3] =
// nvref ;
template <>
const int (*const GenericElement<DataTriangle2>::onWhatBorder)[7] =
    onWhatIsEdge;

template <>
const int (*const GenericElement<DataQuad2>::nvface)[4] = nvfaceQuad;
template <>
const int (*const GenericElement<DataQuad2>::nvedge)[2] = nvedgeQuad;
template <>
const int (*const GenericElement<DataQuad2>::nvhyperFace)[2]       = nvedgeQuad;
template <> const int (*const GenericElement<DataQuad2>::nvadj)[2] = nvedgeQuad;
template <> const int GenericElement<DataQuad2>::nitemdim[4] = {4, 4, 1, 0};
template <>
const int (*const GenericElement<DataQuad2>::commonVertOfEdges)[4] =
    commonVertOfEdgeQuad;
template <>
const int (*const GenericElement<DataQuad2>::onWhatBorder)[9] =
    onWhatIsEdgeQuad;

template <> const int (*const GenericElement<DataSeg2>::nvface)[1] = 0;
template <>
const int (*const GenericElement<DataSeg2>::nvedge)[2] =
    nvedgeSeg; // nvedgeTria ;
// template<> const int (* const GenericElement<DataSeg2>::hyperface)[2] =
// nvedgeSeg; //nvedgeTria ;
template <> const int (*const GenericElement<DataSeg2>::nvadj)[1] = nvadjSeg;
template <> const int GenericElement<DataSeg2>::nitemdim[4] = {2, 1, 0, 0};

// template<> const int (* const GenericElement<DataSeg2>::refElement)[2] =
// nvrefSeg ;

R2 Triangle2::toKref(const R1 &P, int i) const {
   return (1 - P.X()) * R2::KHat[nvedge[i][0]] + P.X() * R2::KHat[nvedge[i][1]];
}
R2 Triangle2::mapToReferenceElement(const R1 &P, int i) const {
   return (1 - P.X()) * R2::KHat[nvedge[i][0]] + P.X() * R2::KHat[nvedge[i][1]];
}

double Triangle2::mesureBord(int i) const {
   return (at(nvedge[i][0]) - at(nvedge[i][1])).norm();
}
