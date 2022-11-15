#include "dataStruct1D.hpp"

// const R1 R1::KHat[2]={R1(0.),R1(1.)};
const std::vector<R1> R1::KHat = {R1(0.), R1(1.)};

// static const int  nvfaceTet[4][3]  = { {2,1,3},{0,2,3},{1,0,3},{0,1,2} };
// static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{0,1},{1,2},{2,3} };

// static const int  nvfaceTria[1][3]  = { {0,1,2} };
// static const int  nvedgeTria[3][2] = { {1,2},{2,0},{0,1}};

static const int nvfaceSeg[1][3] = {{-1, -1, -1}};
static const int nvedgeSeg[1][2] = {{0, 1}};
static const int nvadjSeg[2][1]  = {{0}, {1}};

template <> const int (*const GenericElement<DataPoint1>::nvface)[1] = 0;
template <> const int (*const GenericElement<DataPoint1>::nvedge)[2] = 0;
template <> const int (*const GenericElement<DataPoint1>::nvadj)[1]  = 0;

template <> const int (*const GenericElement<DataSeg1>::nvface)[1] = 0;
template <>
const int (*const GenericElement<DataSeg1>::nvedge)[2] =
    nvedgeSeg; // nvedgeTria ;
template <> const int (*const GenericElement<DataSeg1>::nvadj)[1] = nvadjSeg;

template <> const int GenericElement<DataSeg1>::nitemdim[4] = {2, 1, 0, 0};

static const int onWhatIsVertex[2][3] = {{1, 0, 0},  // sommet  0
                                         {0, 1, 0}}; // sommet 1

template <>
const int (*const GenericElement<DataSeg1>::onWhatBorder)[3] = onWhatIsVertex;
