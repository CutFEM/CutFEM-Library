#include "dataStruct1D.hpp"

const std::vector<R1> R1::KHat = {R1(0.), R1(1.)};

/// @brief Static mumber for Point1
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::nvedge = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::nvface = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::nvhyperFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::edgeOfFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::faceOfEdge = {
    {}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataPoint1>::commonVertOfEdges = {{}};

/// @brief Static mumber for Seg1
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvedge = {{0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvface{{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvhyperFace = {
    {0}, {1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::edgeOfFace = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::faceOfEdge = {{}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataSeg1>::commonVertOfEdges = {{}};
