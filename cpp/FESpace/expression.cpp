/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#include "expression.hpp"
#include "../problem/CutFEM_parameter.hpp"

std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1, double cc) {
    return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst> operator*(double cc, const std::shared_ptr<ExpressionVirtual> &f1) {
    return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                               const Normal_Component_X &cc) {
    return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                               const Normal_Component_Y &cc) {
    return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                               const Normal_Component_Z &cc) {
    return std::make_shared<ExpressionMultConst>(f1, cc);
}

std::shared_ptr<ExpressionAbs> fabs(const std::shared_ptr<ExpressionVirtual> &f1) {
    return std::make_shared<ExpressionAbs>(f1);
}

std::shared_ptr<ExpressionProduct> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                             const std::shared_ptr<ExpressionVirtual> &f2) {
    return std::make_shared<ExpressionProduct>(f1, f2);
}

std::shared_ptr<ExpressionPow> pow(const std::shared_ptr<ExpressionVirtual> &f1, const double nn) {
    return std::make_shared<ExpressionPow>(f1, nn);
}
std::shared_ptr<ExpressionPow> operator^(const std::shared_ptr<ExpressionVirtual> &f1, const double nn) {
    return std::make_shared<ExpressionPow>(f1, nn);
}
std::shared_ptr<ExpressionPow> sqrt(const std::shared_ptr<ExpressionVirtual> &f1) { return pow(f1, 1. / 2); }

std::shared_ptr<ExpressionDivision> operator/(const std::shared_ptr<ExpressionVirtual> &f1,
                                              const std::shared_ptr<ExpressionVirtual> &f2) {
    return std::make_shared<ExpressionDivision>(f1, f2);
}

std::shared_ptr<ExpressionSum> operator+(const std::shared_ptr<ExpressionVirtual> &f1,
                                         const std::shared_ptr<ExpressionVirtual> &f2) {
    return std::make_shared<ExpressionSum>(f1, f2);
}

std::shared_ptr<ExpressionSum> operator-(const std::shared_ptr<ExpressionVirtual> &f1,
                                         const std::shared_ptr<ExpressionVirtual> &f2) {
    return f1 + (-1. * f2);
}

std::shared_ptr<ExpressionNormal2> operator*(const FunFEM<Mesh2> &f1, const Normal &n) {
    return std::make_shared<ExpressionNormal2>(f1, n);
}
std::shared_ptr<ExpressionNormal2> operator*(const FunFEM<Mesh2> &f1, const Tangent &n) {
    return std::make_shared<ExpressionNormal2>(f1, n);
}
std::shared_ptr<ExpressionNormal2> operator*(const FunFEM<Mesh2> &f1, const Conormal &n) {
    return std::make_shared<ExpressionNormal2>(f1, n);
}
std::shared_ptr<ExpressionNormal2Q> operator*(const FunFEM<MeshQuad2> &f1, const Normal &n) {
    return std::make_shared<ExpressionNormal2Q>(f1, n);
}
std::shared_ptr<ExpressionNormal2Q> operator*(const FunFEM<MeshQuad2> &f1, const Tangent &n) {
    return std::make_shared<ExpressionNormal2Q>(f1, n);
}
std::shared_ptr<ExpressionNormal2Q> operator*(const FunFEM<MeshQuad2> &f1, const Conormal &n) {
    return std::make_shared<ExpressionNormal2Q>(f1, n);
}
std::shared_ptr<ExpressionNormal3> operator*(const FunFEM<Mesh3> &f1, const Normal &n) {
    return std::make_shared<ExpressionNormal3>(f1);
}

std::shared_ptr<ExpressionAverage> average(const std::shared_ptr<ExpressionVirtual> &f1, const double kk1,
                                           const double kk2) {
    return std::make_shared<ExpressionAverage>(f1, kk1, kk2);
}
std::shared_ptr<ExpressionAverage> jump(const std::shared_ptr<ExpressionVirtual> &f1, const double kk1,
                                        const double kk2) {
    return std::make_shared<ExpressionAverage>(f1, 1, -1);
}
std::shared_ptr<ExpressionAverage> operator*(double c, const ExpressionAverage &fh) {
    return std::make_shared<ExpressionAverage>(fh.fun1, c * fh.k1, c * fh.k2);
}
std::shared_ptr<ExpressionAverage> operator*(const ExpressionAverage &fh, double c) {
    return std::make_shared<ExpressionAverage>(fh.fun1, c * fh.k1, c * fh.k2);
}


std::vector<std::shared_ptr<ExpressionVirtual>> cross(const Normal &n, const FunFEM<Mesh3> &f1) {

    return {std::make_shared<ExpressionNormalCrossX3>(f1), std::make_shared<ExpressionNormalCrossY3>(f1),
             std::make_shared<ExpressionNormalCrossZ3>(f1)};
}


ExpressionBurgerFlux burgerFlux(const ExpressionVirtual &f1) { return ExpressionBurgerFlux(f1); }
ExpressionNormalBurgerFlux burgerFlux(const ExpressionVirtual &f1, const Normal &n) {
    return ExpressionNormalBurgerFlux(f1);
}


std::shared_ptr<ExpressionDSx3> dxS(const FunFEM<Mesh3> &f1) { return std::make_shared<ExpressionDSx3>(f1); }
std::shared_ptr<ExpressionDSy3> dyS(const FunFEM<Mesh3> &f1) { return std::make_shared<ExpressionDSy3>(f1); }
std::shared_ptr<ExpressionDSz3> dzS(const FunFEM<Mesh3> &f1) { return std::make_shared<ExpressionDSz3>(f1); }
std::shared_ptr<ExpressionDivS3> divS(const FunFEM<Mesh3> &f1) { return std::make_shared<ExpressionDivS3>(f1); }
