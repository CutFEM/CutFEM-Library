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

/**
* @brief A space-time CutFEM for the coupled bulk-surface convection-diffusion equation with a non-linear Langmuir coupling term.
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <array>
#include <iostream>
#ifdef USE_MM_PI
#include "cfmM_PI.hpp"
#endif
#include "finiteElement.hpp"
#include "../num/gnuplot.hpp"
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../num/matlab.hpp"
#include "../num/redirectOutput.hpp"
#include "paraview.hpp"
#include "../problem/AlgoimIntegration.hpp"

using namespace globalVariable; // to access some globally defined constants


namespace Example1 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return ((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17 * 0.17);
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
        return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * (x(0) - 0.5 - 0.28 * sin(M_PI * t)),
                                     2.0 * (x(1) - 0.5 + 0.28 * cos(M_PI * t)));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(pow(2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))), 2) +
                      pow(2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))), 2));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))) / norm,
                  2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))) / norm);
    }
};


// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// Initial solution surface
R fun_uSurfInit(double *P, const int i) {
    double x = P[0], y = P[1];

    return (0.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y) 
    - M_PI / 250 * sin(M_PI * x) * cos(M_PI * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) 
    - M_PI / 250 * cos(M_PI * x) * sin(M_PI * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)))
    / (1.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y));

}


// Exact solution surface
R fun_uSurf(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);

    return (0.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y) * cos(2 * M_PI * t)
    - M_PI / 250 * sin(M_PI * x) * cos(M_PI * y) * cos(2 * M_PI * t) * (x - xc) / sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) 
    - M_PI / 250 * cos(M_PI * x) * sin(M_PI * y) * cos(2 * M_PI * t) * (y - yc) / sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)))
    / (1.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y) * cos(2 * M_PI * t));
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // return M_PI*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*(-4.0/5.0)+((M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI))/1.25E+2-(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*(x-1.0/2.0)*(2.0/5.0)+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*(y-1.0/2.0)*(2.0/5.0);

    return (M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / 125 -
           (4 * M_PI * cos(M_PI * x) * cos(M_PI * y) * sin(2 * M_PI * t)) / 5 -
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * sin(M_PI * y) * (x - 1. / 2)) / 5 +
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * y) * sin(M_PI * x) * (y - 1. / 2)) / 5;
}

// RHS fS surface
R fun_rhsSurf(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (M_PI*1.0/sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.0/pow(cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*4.0+1.5E+1,3.0)*(cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.23046875E+8+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.23046875E+8+cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.26875E+6+cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*2.480625E+7+cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*4.6305E+6+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*7.03125E+7+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.96875E+7+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.96875E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.515625E+7-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.515625E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*3.75E+7+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*5.0E+6+M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.109375E+8-M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.03125E+7-M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.03125E+7+cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*9.375E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*9.375E+6-x*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.0546875E+8-x*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.515625E+8-y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.515625E+8-y*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.0546875E+8-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*9.84375E+7-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.953125E+7+cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.953125E+7+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*9.84375E+7+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.625E+7+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.0E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.625E+7+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.0E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.6875E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.6875E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*2.5E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.5E+6+(x*x)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.0546875E+8+(x*x)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.1640625E+8-(x*x*x)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8+(y*y)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.1640625E+8+(y*y)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.0546875E+8-(y*y*y)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8+pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.480625E+7+pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*8.26875E+6-pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*4.6305E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.05E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.1E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6+M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8-M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.9375E+7-M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*9.375E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.5E+6-M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.5E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*9.375E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.5E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.5E+6+x*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*5.90625E+7+x*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*5.90625E+7+y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.771875E+8+y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*5.90625E+7+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.05E+7-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*9.375E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*9.375E+6-x*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*5.90625E+7-x*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.771875E+8-y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*5.90625E+7-y*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*5.90625E+7+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.65375E+7-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.65375E+7+(x*x)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.53125E+8-(x*x*x)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6875E+8-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.5E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6+(x*x)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.4375E+7+(y*y)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.4375E+7-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.25E+7-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.25E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*2.0E+6-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*2.0E+6+(y*y)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.53125E+8-(y*y*y)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.6875E+8-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.5E+7-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.1025E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*9.84375E+6+x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7-(x*x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8+(y*y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.615E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*6.3E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*5.6E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.9845E+7-pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.7044E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.54E+7-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.12E+6+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.1025E+7+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*5.5125E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.65375E+7-(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.96875E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.96875E+7+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*9.84375E+6+M_PI*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.087E+6+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*7.03125E+7-(x*x)*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.40625E+8-(x*x)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*4.21875E+8+(x*x*x)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*7.03125E+7-(y*y)*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*4.21875E+8-(y*y)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.40625E+8+(y*y*y)*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.9845E+7+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*3.7044E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.54E+7+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*1.12E+6+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*6.615E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*6.3E+6+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*5.6E+5-(x*x)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10-(y*y)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.88E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.25E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.84E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.0E+5+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.5E+7-(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*1.0E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*5.0E+6+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.65375E+7-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*5.5125E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.3075E+7-M_PI*pow(cos(t*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.1025E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*9.84375E+6-(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.96875E+7+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.96875E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),3.0)*sin(y*M_PI)*3.087E+6+M_PI*pow(cos(t*M_PI),3.0)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.174E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7-(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.5E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*5.0E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.0E+7-x*(y*y)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8-(x*x)*y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8-pow(cos(t*M_PI),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.47E+9-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.94E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.88E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*5.88E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.25E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.625E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*5.25E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*7.84E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*1.4E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.8E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*7.0E+5-M_PI*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.1025E+7-M_PI*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.3075E+7+(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*9.84375E+6-M_PI*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.174E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.5E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.5E+6-x*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.65375E+7-(x*x)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*5.90625E+7-y*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*4.96125E+7-(y*y)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.771875E+8-cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.47E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.88E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.94E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.625E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*2.8E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*1.4E+6-x*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*4.96125E+7+(x*x)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.771875E+8-y*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.65375E+7+(y*y)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*5.90625E+7-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*4.6305E+6+pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*4.6305E+6+(x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*6.75E+7-(x*x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*4.5E+7+(x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*6.0E+6-(x*x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6+(x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.25E+7+(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.25E+7+(x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*2.0E+6+(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*2.0E+6+(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*6.75E+7-(y*y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*4.5E+7+(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*6.0E+6-(y*y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*5.5125E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.5435E+6-(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.764E+6+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.568E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.292E+6-pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*9.8784E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.704E+5-pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*8.7808E+4+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*5.5125E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.65375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*1.5435E+6-(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.087E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.109375E+8-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.8125E+8+(x*x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.40625E+8-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.109375E+8+(y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.8125E+8-(y*y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.40625E+8-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.875E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.875E+7+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*5.292E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*9.8784E+5+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*4.704E+5+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*8.7808E+4+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.764E+6+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.568E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.94E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*8.232E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*3.92E+5+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0976E+5-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.8125E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.8125E+7+(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*1.0E+7-(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.65375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*3.087E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.8125E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.8125E+7+(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.0E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.205E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.35E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*2.94E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.116E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*8.232E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*7.84E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*3.92E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*1.0976E+5+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.34375E+9-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.34375E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.4375E+8-x*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.4375E+7-y*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.4375E+7-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.4375E+8-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.35E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.205E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*4.116E+5-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),3.0)*7.84E+5-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.3625E+7-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*6.825E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.96875E+7+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7-x*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*4.21875E+8+x*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.40625E+8+x*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*4.21875E+8+y*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8+y*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.40625E+8+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*6.825E+7+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.3625E+7-x*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10-y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+y*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.5E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8+M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.9375E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.96875E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.5E+7+x*y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8+x*y*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-x*(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*4.5E+7-x*(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6-(x*x)*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*4.5E+7-(x*x)*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.1025E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*5.90625E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.087E+6-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.1025E+7+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.96875E+7+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+(x*x)*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*8.82E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.6464E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*5.88E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*6.25E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*6.25E+8-x*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*3.528E+6-x*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*3.136E+5-(x*x)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.26E+7-(x*x)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.12E+6-y*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0584E+7-(y*y)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*3.78E+7-y*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*9.408E+5-(y*y)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*3.36E+6-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.1025E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.1025E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*1.96875E+7-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.1025E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.3075E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.18125E+8+(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*5.90625E+7-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*3.087E+6-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.75E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.75E+7-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.0584E+7+(x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.78E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*9.408E+5+(x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*3.36E+6+(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.26E+7+(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*1.12E+6-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.528E+6-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.136E+5+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.88E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*3.15E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.6464E+6-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.84E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*4.2E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.1952E+5-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.5435E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*1.5435E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.88E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.05E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.84E+5+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.4E+6-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.3075E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.18125E+8+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.625E+7-x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.625E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.625E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.625E+7-x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7+x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7-(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7+(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.0E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*2.0E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.1025E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.771875E+8-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.1E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.05E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*9.8784E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*8.7808E+4-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*9.8784E+5-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.7808E+4+x*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+y*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.41E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.47E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*5.88E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.05E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*8.232E+5-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.05E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.568E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*7.84E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*5.6E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.4E+6+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.087E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.087E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.41E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.47E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*5.88E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.05E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*3.15E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*1.6464E+6-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*7.84E+5+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*4.2E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*2.1952E+5-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.771875E+8-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.5E+6-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.05E+7-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.5E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.1E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+6+(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6875E+8+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*8.232E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*8.232E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*1.0976E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.0976E+5+x*y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.6875E+8-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.47E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.41E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.05E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.47E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.41E+6+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*8.232E+5-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.05E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),3.0)*1.568E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*5.6E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.96E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*7.84E+5+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*7.875E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+5+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0976E+5-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+x*y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*4.725E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.725E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*4.725E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.4175E+8+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.116E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.116E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),3.0)*4.3904E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*4.3904E+5+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*7.875E+7+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.88E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0976E+5-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.18125E+8-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.96E+5-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*7.84E+5-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+x*y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*8.4375E+8-x*y*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8-x*y*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.4175E+8-x*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*4.725E+7-y*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*4.725E+7-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*4.725E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*4.2E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*5.6E+6+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.3625E+8-x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7-x*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.18125E+8+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-y*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.3625E+8-y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7+M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.323E+7-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.323E+7-y*cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.05E+10-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*4.2E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*5.6E+6-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.1025E+7+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.1025E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.1E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.2E+7+x*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+x*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.3625E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-y*M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.3625E+8+y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+y*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.75E+7+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.625E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*5.625E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.125E+8-(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.125E+8-(x*x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7+(y*y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.0E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.0E+7+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.75E+7+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.0E+6-x*y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.18125E+8+x*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.05E+10-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.1E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*6.615E+7+M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.205E+7+M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.205E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.82E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.575E+7+M_PI*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6464E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.6E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.875E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.94E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.4E+6+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*7.5E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+7+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.875E+7-(x*x)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+x*y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.18125E+8-(y*y)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*2.94E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7-M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.8E+6-M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.4E+6+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.82E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*5.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.575E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(y*M_PI)*1.6464E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.40625E+10-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-(x*x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.40625E+10+(y*y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.3075E+7-pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8+y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.3075E+7-x*(y*y)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6875E+8+x*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*4.5E+7+x*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6-(x*x)*y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.6875E+8+x*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*4.5E+7+x*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.96875E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8-x*(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+(x*x)*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*5.25E+6+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*2.1E+7-cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8-x*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.323E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.26E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.12E+6-(x*x)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.725E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.26E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.26E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.12E+6+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.12E+6-y*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.969E+7-(y*y)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.4175E+8+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*3.78E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*3.36E+6-x*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.1025E+7-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*1.96875E+7+(x*x)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.3075E+7+y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.18125E+8-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7+(y*y)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.18125E+8+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+x*(y*y)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8+(x*x)*y*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.969E+7+(x*x)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.4175E+8-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.78E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*3.36E+6+(y*y)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*4.725E+7-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.26E+7-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.26E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*1.12E+6-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*1.12E+6-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.323E+7-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.26E+7-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*1.12E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.05E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.4E+6-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.3075E+7+x*M_PI*pow(cos(t*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.205E+7+x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7+x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.18125E+8+(x*x)*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.18125E+8-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.625E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.625E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*2.0E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.0E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.1025E+7+y*M_PI*pow(cos(t*M_PI),2.0)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.615E+7-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*7.875E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7+(y*y)*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.3625E+8+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.7044E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.528E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*3.136E+5-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.7044E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*3.528E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*3.136E+5-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7+x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.25E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.05E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*5.6E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*5.6E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.4E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.087E+6+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.087E+6+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.1025E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.1025E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.05E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*5.6E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6+x*M_PI*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.615E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*7.875E+7-(x*x)*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.3625E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8+y*M_PI*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.205E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7-(y*y)*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*5.625E+7-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.75E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.0E+6-(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7-(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.0E+6-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.875E+7-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.625E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*5.625E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.0E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.5E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.0E+7-(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.5E+7-(x*x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+(y*y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7-M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.875E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.625E+7-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.75E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.0E+6-(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.0E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.625E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.05E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*5.6E+6+M_PI*cos(t*M_PI)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.174E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.174E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.25E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*5.6E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*5.6E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*4.41E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.568E+6-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.176E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.0E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6+M_PI*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.1952E+5+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.125E+8-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+7+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.41E+6-(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.232E+5-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*7.84E+5+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*3.92E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.4E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.4E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*7.5E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.125E+8-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*7.5E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.47E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.47E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*1.568E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*1.568E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*4.41E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*8.232E+5-M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.84E+5-M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.92E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.4E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.4E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*8.82E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*4.41E+6-(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.6464E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.568E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.176E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*7.0E+5-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(y*M_PI)*2.1952E+5-x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*4.41E+7-y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*4.41E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.15E+7+(x*x)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.1E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6-x*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.88E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7+y*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.764E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.15E+7+(y*y)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.3E+7+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.8E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.12E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*8.4E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+8+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.8E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.5E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.6E+6+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*2.8E+6-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.764E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.15E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*5.88E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*6.3E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.12E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*8.4E+6+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*2.8E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*6.3E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.15E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.1E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6-x*(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+(x*x)*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.625E+7-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*7.5E+7+(x*x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.75E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.625E+7+(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*7.5E+7-(y*y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.75E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.6464E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*8.82E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*7.84E+5-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*6.3E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.6464E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*8.82E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*7.84E+5-x*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+y*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+(y*y)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.47E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5-x*y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.52E+7-x*y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*2.24E+6+x*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.875E+7+x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.176E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.47E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5+x*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.52E+7+x*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*2.24E+6-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6+y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7-x*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*7.875E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*7.056E+6+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*7.056E+6+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*6.272E+5+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*6.272E+5-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.205E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.15E+7-x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.1E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.12E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*2.205E+7-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7+x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.125E+8+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.5E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.0E+7+(x*x)*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+x*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.75E+7+y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.75E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.125E+8-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.0E+7-x*(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.125E+8-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.125E+8-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*3.15E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*1.12E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.82E+6+x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.4E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.6E+6+(x*x)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.05E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7-x*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*7.84E+5-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.6E+6+y*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.352E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.8E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.6E+6+(y*y)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*8.4E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+8-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.82E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.5E+8-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*3.136E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.88E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*3.136E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*8.82E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*2.352E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*5.6E+6+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*7.84E+5-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*5.6E+6+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*8.4E+6-M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.05E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.764E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.82E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*6.3E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.4E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*5.6E+6+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*2.205E+7-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.4696E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*2.1952E+5-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*3.92E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.568E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.764E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*6.3E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*2.4696E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.1952E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*3.92E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.568E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.176E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*2.205E+7-(x*x)*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*3.2928E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*3.2928E+6-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8+x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*3.087E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.176E+7-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.568E+6-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.087E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7+(x*x)*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.1E+7+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.1E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.3E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.1E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*8.82E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.94E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*1.176E+7-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.15E+7+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.1E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*1.568E+6-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.1E+7-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.125E+8+(x*x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7-x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.5E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.0E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.5E+7+x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7+(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*6.3E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.1E+7-(x*x)*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.75E+7-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.75E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.0E+7-x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7-(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.1E+7-(y*y)*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.125E+8+(y*y*y)*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7+x*y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.875E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.875E+7-y*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*1.6464E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*2.1952E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.6464E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.1952E+5-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.94E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*8.82E+6+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*3.15E+7-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*3.92E+5+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.568E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*5.88E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.94E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*3.92E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.6E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.8E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.1952E+5-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*5.88E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*8.82E+6+M_PI*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.6464E+6+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+x*y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.625E+6+x*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-x*y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*9.45E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*8.232E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*8.232E+5-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.88E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.176E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6-x*y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.875E+7-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*8.82E+6-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.6464E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.92E+5-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.176E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*5.6E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.1952E+5-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.8E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.94E+6+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.92E+5+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.568E+6+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.625E+6+x*y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*9.45E+7+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.0976E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*4.3904E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.0976E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*4.3904E+5-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.575E+8-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*7.875E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-y*M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.75E+9-(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.75E+9-(x*x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+(y*y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.646E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.646E+7-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*2.205E+7+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.205E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-x*y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.575E+8-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.47E+9-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.125E+8+(x*x)*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7-x*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.125E+8-x*(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7+x*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.764E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.2E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*4.2E+7+x*(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7+(x*x)*y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.176E+7+x*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*5.88E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.84E+5+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.6E+6+(x*x)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.352E+6-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*8.4E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.176E+7+y*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.764E+7+(y*y)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.3E+7-x*y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.1E+7+x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.764E+7-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.3E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.352E+6-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*8.4E+6-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.84E+5+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6+y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*5.88E+6-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.1E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.6464E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.1952E+5-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.1952E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.568E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.568E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.6464E+6+(x*x)*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+(y*y)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.3E+7-x*y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.6E+6-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*8.4E+7-x*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*5.6E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*6.3E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.568E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.764E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.764E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.568E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*8.4E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6-(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.725E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*2.352E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*2.352E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.725E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7-x*y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7+x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.25E+8+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.6E+6+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.12E+7-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*8.4E+6-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.3E+7+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*6.3E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.3E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*8.4E+6+x*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7+y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.12E+7+y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*6.3E+7+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*5.88E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.568E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*3.136E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.568E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.136E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*5.88E+6+(x*x)*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*4.2E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-x*y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*4.2E+7+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.176E+7-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.176E+7-x*y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7))/(x*-1.5625E+7-y*1.5625E+7-cos(t*M_PI)*4.375E+6+sin(t*M_PI)*4.375E+6+y*cos(t*M_PI)*8.75E+6-x*sin(t*M_PI)*8.75E+6+(x*x)*1.5625E+7+(y*y)*1.5625E+7+pow(cos(t*M_PI),2.0)*1.225E+6+pow(sin(t*M_PI),2.0)*1.225E+6+7.8125E+6);;
}

R fun_one(double *P, const int i) { return 1.; }

} // namespace Example1


namespace Kex {

    template <int N> struct Levelset {

        double t;

        // level set function
        template <typename V> typename V::value_type operator()(const V &P) const {
            R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
            return -(((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17));
        }

        // gradient of level set function
        template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

            return algoim::uvector<T, N>(-2.0 * (x(0) - 0.5 - 0.28 * sin(M_PI * t)),
                                        -2.0 * (x(1) - 0.5 + 0.28 * cos(M_PI * t)));
        }

        // normal = grad(phi)/norm(grad(phi))
        R2 normal(std::span<double> P) const {
            R norm = sqrt(pow(2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))), 2) +
                        pow(2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))), 2));
            // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
            return R2(-2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))) / norm,
                    -2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))) / norm);
        }
    };

    R fun_one(double *P, const int i) { return 1.; }

    // Level-set function
    double fun_levelSet(double *P, const int i, const R t) {
        R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
        return -(sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17 - Epsilon);
    }

    // Level-set function initial
    double fun_levelSet(double *P, const int i) {
        return -(sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17 - Epsilon);
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * cos(pi * y) * sin(pi * x) * (2 * cos(pi * t) * cos(pi * t) - 1) * ((7 * sin(pi * t)) / 25 - x + 0.5)) /
                (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2))) -
            (pi * cos(pi * x) * sin(pi * y) * (2 * cos(pi * t) * cos(pi * t) - 1) * (y + (7 * cos(pi * t)) / 25 - 0.5)) /
                (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2)));
    }

    // Velocity field
    R fun_velocity(double *P, const int i) {
        if (i == 0)
            return M_PI * (0.5 - P[1]);
        else
            return M_PI * (P[0] - 0.5);
    }

    // Initial solution bulk
    R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

    // Exact solution bulk
    R fun_uBulk(double *P, const int i, const R t) {
        return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
    }

    R fun_uBulkD(double *P, const int i, const int d, const R t) {
        return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
    }

    // Initial solution surface
    R fun_uSurfInit(double *P, const int i) {
        double x = P[0], y = P[1];

        return (0.5 + 0.4 * cos(pi * x) * cos(pi * y) +
            pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) +
            pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)))
            /(1.5 + 0.4 * cos(pi * x) * cos(pi * y));
    }

    // Exact solution surface
    R fun_uSurf(double *P, const int i, const R t) {
        double x = P[0], y = P[1];

        R xc = 0.5 + 0.28 * sin(pi * t), yc = 0.5 - 0.28 * cos(pi * t);

        return (0.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t) +
            pi / 250 * sin(pi * x) * cos(pi * y) * cos(2 * pi * t) * (x - xc) /
                sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) +
            pi / 250 * cos(pi * x) * sin(pi * y) * cos(2 * pi * t) * (y - yc) /
                sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)))
                /(1.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t));
    }

    

    // RHS fB bulk
    R fun_rhsBulk(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * pi * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 125 -
            (4 * pi * cos(pi * x) * cos(pi * y) * sin(2 * pi * t)) / 5 -
            (2 * pi * pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (x - 1. / 2)) / 5 +
            (2 * pi * pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (y - 1. / 2)) / 5;
    }



    // RHS fS surface
    R fun_rhsSurf(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        // return -(pi*(17578125*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 17578125*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 8268750*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 24806250*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 4630500*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)+ 70312500*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 4687500*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 4687500*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 35156250*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 35156250*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 37500000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 5000000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 210937500*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 70312500*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 70312500*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 4687500000*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 4687500000*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 9375000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 105468750*x*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 70312500*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 70312500*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 105468750*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 19687500*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 29531250*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 29531250*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 19687500*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 6250000*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1000000*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 6250000*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 1000000*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 4687500*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 4687500*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 105468750*x*x*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 316406250*x*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*pow(x,3)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 316406250*y*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 105468750*y*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*pow(y,3)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 24806250*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 8268750*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 4630500*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 10500000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 21000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 21000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 2800000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 118125000*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 39375000*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 2500000*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 2500000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 2500000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 2500000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 59062500*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 59062500*x*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 177187500*y*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 59062500*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 2625000000*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 21000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 10500000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 59062500*x*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 177187500*x*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 59062500*y*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 59062500*y*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 140625000*x*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 140625000*y*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 16537500*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 16537500*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 253125000*x*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 168750000*pow(x,3)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 35000000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 4000000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 84375000*x*x*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 84375000*y*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 22500000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 22500000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2000000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 2000000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 253125000*y*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 168750000*pow(y,3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 35000000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 4000000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 11025000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 9843750*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 70312500*x*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*x*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*pow(x,3)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 70312500*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 70312500*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 210937500*y*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 140625000*pow(y,3)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 6615000*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 6300000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 560000*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 19845000*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 3704400*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 9800000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 1120000*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 75000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 10000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 11025000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 5512500*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 16537500*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 19687500*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 19687500*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 9843750*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 3087000*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 75000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 10000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 70312500*x*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 140625000*x*x*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 421875000*x*x*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 281250000*pow(x,3)*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 70312500*y*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 421875000*y*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 140625000*y*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 281250000*pow(y,3)*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 19845000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 3704400*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 9800000*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 1120000*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi) 
        // - 6615000*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 6300000*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 560000*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi)+ 18750000000*x*x*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 18750000000*y*y*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5880000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 5250000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 784000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 700000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 18750000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 75000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 37500000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 10000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 5000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 16537500*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 5512500*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 33075000*pi*pow(cos(t*pi),2)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 11025000*pi*pow(cos(t*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 9843750*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 19687500*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 19687500*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 3087000*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),3)*sin(y*pi)+ 6174000*pi*pow(cos(t*pi),3)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 18750000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 37500000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 75000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 5000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 10000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 210937500*x*y*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*x*x*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 1470000000*pow(cos(t*pi),2)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2940000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 5880000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 5880000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 5250000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 2625000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 5250000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 784000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 1400000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 2800000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 700000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 11025000*pi*cos(x*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 33075000*pi*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi)+ 9843750*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 6174000*pi*cos(y*pi)*pow(sin(t*pi),3)*sin(2*t*pi)*sin(x*pi)+ 1250000*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1250000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 16537500*x*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 59062500*x*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 49612500*y*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 177187500*y*y*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 1470000000*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5880000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 2940000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 2625000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 5250000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 2800000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3) 
        // - 1400000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 49612500*x*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 177187500*x*x*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 16537500*y*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 59062500*y*y*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 4630500*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 4630500*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 67500000*x*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 45000000*pow(x,3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 6000000*x*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 4000000*pow(x,3)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 22500000*x*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 22500000*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 2000000*x*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 2000000*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 67500000*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 45000000*pow(y,3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 6000000*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 4000000*pow(y,3)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 5512500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 1543500*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 70312500*x*x*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 210937500*x*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*pow(x,3)*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*y*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*y*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*pow(y,3)*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 1764000*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 156800*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 5292000*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 987840*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 470400*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 87808*pow(cos(t*pi),3)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 5512500*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 16537500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 5512500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 1543500*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3) 
        // - 3087000*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 210937500*x*x*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 281250000*pow(x,3)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 140625000*pow(x,4)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 210937500*y*y*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 281250000*pow(y,3)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 140625000*pow(y,4)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 18750000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 18750000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 5292000*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 987840*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 470400*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 87808*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 1764000*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 156800*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 2940000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 823200*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 392000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 109760*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 28125000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 28125000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 18750000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 10000000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 5512500*pi*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 16537500*pi*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 5512500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 3087000*pi*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)+ 1543500*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 28125000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 28125000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 18750000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 10000000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 4687500000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1250000000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2205000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 735000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 2940000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 411600*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 823200*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3) 
        // - 784000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 392000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 109760*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),3) 
        // - 1250000000*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5512500*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 1543500*pi*pi*cos(2*t*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 140625000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 140625000*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 2343750000*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2343750000*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 93750000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 84375000*x*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 84375000*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 93750000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 735000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 2205000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 411600*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*pow(sin(x*pi),2) 
        // - 784000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(t*pi),2)*pow(sin(y*pi),3) 
        // - 39375000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 140625000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 70312500*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 1250000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 70312500*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 140625000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 4687500000*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 23625000*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 26250000*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 75000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 10000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 19687500*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 39375000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 39375000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 75000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 10000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 421875000*x*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 140625000*x*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 281250000*x*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 1250000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1250000000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 421875000*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 281250000*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 140625000*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 26250000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 23625000*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 9375000000*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000000*x*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000000*y*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 2800000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 75000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 37500000*x*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 118125000*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 78750000*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 39375000*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 39375000*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 19687500*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 37500000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 75000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 210937500*x*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 210937500*x*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 2625000000*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5250000000*cos(t*pi)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 45000000*x*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 4000000*x*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 45000000*x*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 4000000*x*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 11025000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 59062500*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 3087000*x*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 39375000*pow(x,3)*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 11025000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 19687500*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 140625000*x*y*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 140625000*x*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 8820000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(sin(t*pi),2)*sin(x*pi)*pow(sin(y*pi),2)+ 1646400*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(sin(t*pi),3)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 5880000*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(x*pi),2)*sin(y*pi)+ 625000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 625000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 3528000*x*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 313600*x*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 12600000*x*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 1120000*x*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 10584000*y*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 37800000*y*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 940800*y*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 3360000*y*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 11025000*x*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 11025000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 39375000*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 19687500*x*x*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 11025000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2) 
        // - 33075000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 118125000*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 59062500*y*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 3087000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3) 
        // - 39375000*pow(y,3)*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 37500000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 75000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 75000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 37500000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 10584000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 37800000*x*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 940800*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 3360000*x*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi) 
        // - 12600000*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 1120000*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 3528000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 313600*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(y*pi)+ 5880000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 31500000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 1646400*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 21000000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 784000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 4200000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 219520*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 2800000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 1543500*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 1543500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 5880000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 10500000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 784000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 1400000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 33075000*x*pi*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 11025000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 19687500*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 118125000*x*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 56250000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 56250000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 56250000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 56250000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 37500000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 37500000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 37500000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 37500000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 20000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 20000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 11025000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 33075000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 177187500*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 39375000*y*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 3087000*y*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 118125000*pow(y,3)*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 21000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 10500000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 9375000000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 987840*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 87808*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 987840*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 87808*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 5000000000*x*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5000000000*y*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 4410000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 1470000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 5880000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 15750000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 10500000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 823200*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 10500000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 1568000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 784000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 5600000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 1400000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 3087000*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 3087000*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 4410000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 1470000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 5880000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 15750000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 10500000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 31500000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 1646400*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3) 
        // - 21000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 784000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 4200000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 219520*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),3) 
        // - 2800000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 2500000000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 33075000*x*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 177187500*x*x*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 3087000*x*pi*pi*cos(2*t*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 118125000*pow(x,3)*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 11025000*y*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 19687500*y*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 7500000*x*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 10500000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 7500000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 10000000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 7500000*y*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 10000000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 5000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 2625000000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 21000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 9375000000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 4687500000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 7500000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 4687500000*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 9375000000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1400000000*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 168750000*x*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 823200*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 823200*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 109760*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 109760*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 168750000*x*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 700000000*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1470000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 4410000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 10500000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 15750000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 1543500*pi*pi*cos(t*pi)*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 1543500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 1470000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 4410000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 15750000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 823200*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*pow(sin(x*pi),2) 
        // - 10500000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(t*pi),2)*pow(sin(y*pi),3) 
        // - 5600000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 196000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 784000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 78750000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 588000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 109760*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 2500000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 140625000*x*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 140625000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 1312500000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1312500000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2625000000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*x*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 47250000*x*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 47250000*x*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 47250000*y*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 141750000*y*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 700000000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 411600*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 411600*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 439040*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(t*pi),2)*pow(sin(y*pi),3) 
        // - 439040*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 39375000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 78750000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 588000*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 109760*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),3)*pow(sin(x*pi),2)+ 700000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 39375000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 118125000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 196000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 784000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 5000000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 843750000*x*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 281250000*x*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 281250000*x*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 2625000000*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1312500000*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1312500000*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5000000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 141750000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 47250000*x*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 47250000*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 47250000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 42000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 5600000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 236250000*x*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 78750000*x*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*x*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 118125000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 39375000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 1400000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 700000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 236250000*y*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*y*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 39375000*y*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 39375000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 2625000000*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 13230000*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 13230000*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 10500000000*y*cos(t*pi)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 5600000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 11025000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 11025000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 21000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 42000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 78750000*x*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 236250000*x*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 1400000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 236250000*y*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 78750000*y*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 78750000*y*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 37500000*x*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 10000000*x*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 5000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 56250000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 56250000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 112500000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 112500000*y*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 75000000*pow(x,3)*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 75000000*pow(y,3)*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 20000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 10000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 20000000*y*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 37500000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 10000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 118125000*x*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 10500000000*x*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 21000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 11025000*x*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 66150000*pi*cos(t*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 22050000*pi*cos(t*pi)*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 22050000*pi*cos(t*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 33075000*y*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 8820000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 15750000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 1646400*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 5600000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 2800000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 18750000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 2940000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 2800000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 1400000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 9375000000*x*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 75000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 75000000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 9375000000*y*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 5000000000*x*x*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 118125000*x*y*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 5000000000*y*y*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 33075000*x*x*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 11025000*y*y*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 2940000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 21000000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 2800000*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1400000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 735000000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 8820000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 5250000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 15750000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 1646400*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),3)*sin(y*pi)+ 5600000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 2800000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 14062500000*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 4687500000*x*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*pow(x,3)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 4687500000*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 14062500000*y*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 9375000000*pow(y,3)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 33075000*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 392000000*pow(cos(t*pi),2)*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 33075000*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 168750000*x*y*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 45000000*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 4000000*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 168750000*x*x*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 45000000*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 4000000*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 39375000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 19687500*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 2500000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 140625000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 140625000*x*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*x*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 140625000*x*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 5250000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 735000000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 367500000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 367500000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 392000000*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 13230000*x*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 12600000*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1120000*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 47250000*x*x*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 12600000*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 12600000*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1120000*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 1120000*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 39690000*y*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 141750000*y*y*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 37800000*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 3360000*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 11025000*x*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 19687500*x*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 39375000*x*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 196000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 33075000*y*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 118125000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 39375000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 118125000*y*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 5000000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 281250000*x*y*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 281250000*x*x*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 367500000*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 367500000*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5000000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 39690000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 141750000*x*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 37800000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 3360000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi) 
        // - 47250000*y*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 12600000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 12600000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 1120000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi)+ 1120000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 13230000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 12600000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 1120000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 21000000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 2800000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 10500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 1400000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 33075000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 22050000*x*pi*pow(cos(t*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 19687500*x*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 39375000*x*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 118125000*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 78750000*x*x*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 118125000*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 56250000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 56250000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 392000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 196000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 11025000*y*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)+ 66150000*y*pi*pow(cos(t*pi),2)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 39375000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 39375000*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 236250000*y*y*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 39375000*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 3704400*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 3528000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 313600*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 3704400*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 3528000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 313600*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 15750000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 5250000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 10500000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 1400000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 3087000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)+ 3087000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 11025000*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 11025000*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 15750000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 10500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 21000000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 5600000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 2800000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 66150000*x*pi*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 78750000*x*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 236250000*x*x*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 392000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 22050000*y*pi*cos(x*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi)+ 19687500*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 78750000*y*y*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 56250000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 37500000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 5000000*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 5000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 10000000*x*x*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 18750000*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 56250000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 56250000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 5000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 10000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 5000000*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 20000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 15000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 20000000*y*y*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 15000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 10000000*pow(x,3)*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 10000000*pow(y,3)*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 18750000*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 56250000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 37500000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 5000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 10000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 56250000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 10500000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 15750000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 5600000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 6174000*pi*cos(t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 6174000*pi*pow(cos(t*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 5250000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 15750000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3) 
        // - 4410000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 1568000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 1176000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 700000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 2800000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 219520*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 112500000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 75000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 4410000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 823200*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 784000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 392000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 1400000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 1400000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 75000000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 75000000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 112500000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 75000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 1470000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 1470000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3) 
        // - 5880000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 4410000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 823200*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 784000*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 392000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi)+ 1400000*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1400000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 8820000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 4410000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 1646400*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 1568000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 1176000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 700000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 2800000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 219520*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(y*pi) 
        // - 44100000*x*pi*cos(t*pi)*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 44100000*y*pi*cos(t*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 31500000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 21000000*x*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 2800000*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 5880000*x*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 31500000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 17640000*y*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 31500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 63000000*y*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 2800000*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 11200000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 8400000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 150000000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 31500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 5600000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 2800000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 98000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 150000000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 98000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 31500000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 5600000*x*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 2800000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 17640000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 42000000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 31500000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 5880000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 42000000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 31500000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 63000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 11200000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 8400000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 2800000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 63000000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 31500000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 21000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 2800000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 9375000000*x*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 9375000000*x*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 56250000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 75000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 37500000*pow(x,4)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 56250000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 75000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 37500000*pow(y,4)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 39375000*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 1646400*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 8820000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 784000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 63000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 1646400*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 8820000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 784000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 735000000*x*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 42000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 42000000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 735000000*y*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5250000000*y*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 42000000*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 1470000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 25200000*x*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2240000*x*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 78750000*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 39375000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 11760000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 11760000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 735000000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5250000000*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 735000000*y*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1470000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 25200000*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 2240000*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 22050000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 21000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 2800000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 22050000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 39375000*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 78750000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 7056000*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 7056000*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 627200*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 627200*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 1250000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 22050000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 31500000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 21000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 21000000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 11200000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 2800000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 22050000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 39375000*x*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 75000000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 112500000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 15000000*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 10000000*x*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 37500000*x*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 37500000*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 112500000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 15000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 10000000*x*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 75000000*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 112500000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 112500000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 11760000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 1568000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 11760000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 1568000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 21000000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 31500000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 11200000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 8820000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 1400000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2800000*x*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 10500000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 31500000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 784000*x*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 4200000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 2352000*y*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2800000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5600000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 8400000*y*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 21000000*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 150000000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 8820000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 31500000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 4200000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 150000000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 31500000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 5880000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 3136000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 5880000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 3136000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 11760000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 8820000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 31500000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 4200000*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 21000000*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 42000000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 31500000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 2352000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 2800000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 5600000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 784000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 4200000*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 8400000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 10500000*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 17640000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 8820000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 63000000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 1400000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 2800000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 31500000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 22050000*x*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 39375000*x*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 2469600*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 219520*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 392000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 17640000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(sin(t*pi),2)*sin(x*pi)*pow(sin(y*pi),2)+ 63000000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 2469600*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 219520*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 392000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 42000000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 11760000*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 1250000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1250000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 22050000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2) 
        // - 39375000*x*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 3292800*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 3292800*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 75000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 350000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 75000000*x*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 350000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 3087000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 22050000*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 11760000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 21000000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 1568000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 2800000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 3087000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 22050000*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 39375000*x*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 21000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 21000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 63000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 21000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 350000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 350000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 8820000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 2940000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 11760000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2) 
        // - 31500000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 21000000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 21000000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 1568000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2) 
        // - 2800000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 39375000*x*y*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 21000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 112500000*x*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 75000000*pow(x,3)*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 15000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 20000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 15000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 63000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 21000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 37500000*x*x*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 37500000*y*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 15000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 15000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 20000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 10000000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 10000000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5250000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 112500000*y*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 75000000*pow(y,3)*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 9375000000*x*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 18750000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 2800000000*y*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1646400*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 11760000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 219520*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 1568000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 1646400*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 11760000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 219520*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 1568000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 3087000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 2940000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 8820000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 21000000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 31500000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 3087000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 392000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 1568000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 5880000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 2940000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 1176000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 392000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 4200000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5600000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 2800000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 219520*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 2800000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5880000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 8820000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 1646400*pi*pow(cos(t*pi),3)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 2625000000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1176000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 4200000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 2625000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000000*x*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2625000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 2800000000*x*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 94500000*x*y*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 823200*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 5880000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 823200*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 5880000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1176000*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 4200000*x*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 78750000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 8820000*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 1646400*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(2*t*pi)*sin(x*pi)+ 392000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 1176000*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 2800000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 4200000*y*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 5600000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 219520*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),3)*pow(sin(x*pi),2) 
        // - 2800000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1400000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2940000*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 5250000000*x*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2625000000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 392000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 1568000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 2625000000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2625000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 94500000*x*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 109760*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 439040*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 109760*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 439040*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 1400000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 157500000*x*y*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 78750000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 2800000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*y*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1250000000*x*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 3750000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 3750000000*y*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*pow(x,3)*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*pow(y,3)*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1250000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 26460000*x*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 26460000*y*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 22050000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 22050000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 2800000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 157500000*x*y*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 1470000000*pi*cos(t*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 112500000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 75000000*x*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 10000000*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 98000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 112500000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 75000000*x*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 10000000*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 98000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2940000*x*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 17640000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 8820000*y*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 2500000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*x*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 8820000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 700000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2940000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 700000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 700000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 700000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 75000000*x*y*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 75000000*x*x*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 11760000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 5880000*x*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 784000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 5600000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 21000000*x*x*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 2352000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 8400000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 11760000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 17640000*y*pi*pow(cos(t*pi),2)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 63000000*y*y*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 5250000000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 21000000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 17640000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 63000000*x*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 2352000*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 8400000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 21000000*y*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 784000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 5600000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 5880000*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi)+ 5250000000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 5250000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 1646400*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 219520*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 219520*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 1568000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 1646400*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 2500000000*x*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*x*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 196000000*x*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 196000000*y*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1400000000*y*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 63000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 5600000*x*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 196000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 196000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1400000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 84000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 5600000*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 63000000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 1568000*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 17640000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 17640000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 1568000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 84000000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 2940000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 5250000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 8820000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 47250000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*y*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 31500000*pow(y,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 23520000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 23520000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 8820000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 47250000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 31500000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 2940000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 5250000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 2500000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 75000000*x*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 75000000*x*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 225000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 700000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 5600000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 11200000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 21000000*x*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 21000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 8400000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 63000000*y*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 63000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 700000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 63000000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 8400000*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 21000000*x*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 21000000*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 5600000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 11200000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 21000000*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 63000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 5880000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 1568000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 3136000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 1568000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 3136000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 5880000*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 10500000*x*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 10500000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 5000000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 42000000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 1400000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 392000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 11760000*x*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 11760000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 1400000000*x*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 - x+0.5),2))) 
        // - 10500000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 10500000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)))
        // /(12500*sqrt(pow((y+ (7*cos(pi*t))/25 -0.5),2)+ pow(((7*sin(pi*t))/25 - x+0.5),2))
        // *pow((4*cos(2*pi*t)*cos(pi*x)*cos(pi*y)+ 15),3)*(350*sin(t*pi) - 1250*y- 350*cos(t*pi) 
        // - 1250*x+ 700*y*cos(t*pi) - 700*x*sin(t*pi)+ 1250*x*x+ 1250*y*y+ 98
        // *pow(cos(t*pi),2)+ 98*pow(sin(t*pi),2)+ 625));


        return -(pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                 1) *
                   ((pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) /
                         (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                     1) *
                        (-((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                               (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) +
                           (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) -
                           (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                            pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                           2.5)) +
                           (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                            pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                           2.5)) +
                           (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                            (2 * y + (14 * cos(t * pi)) / 25 - 1) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) +
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                            (2 * y + (14 * cos(t * pi)) / 25 - 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5))) /
                             ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                         (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                         (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(x * pi) * sin(y * pi) *
                          sin(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                         (4 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                          (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) -
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)))) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) +
                    ((2 * y + (14 * cos(t * pi)) / 25 - 1) /
                         (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                     ((2 * y + (14 * cos(t * pi)) / 25 - 1) * pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2)) /
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             2)) *
                        ((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                          (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                           (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                           (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) -
                          (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                           (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5))) /
                             ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                         (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) +
                    ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                    ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)))) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                    ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (2 * y + (14 * cos(t * pi)) / 25 - 1) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                            2)) -
               ((4 * pi * cos(x * pi) * cos(y * pi) * sin(2 * t * pi)) / 5 -
                (pi * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                    (125 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                (7 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                    (6250 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                (7 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi)) /
                    (6250 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                (pi * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                    (125 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 ((14 * pi * sin(t * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) / 25 -
                  (14 * pi * cos(t * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) / 25) *
                 (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                    (500 *
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         1.5)) +
                (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 ((14 * pi * sin(t * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) / 25 -
                  (14 * pi * cos(t * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) / 25) *
                 (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                    (500 *
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         1.5))) /
                   ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
               (pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) /
                    (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                1) *
                   ((pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) /
                         (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                     1) *
                        (((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                          (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                              (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                          (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                           (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) +
                          (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                           (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                          (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) -
                          (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) -
                          (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                           pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          2.5)) +
                          (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                           (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) +
                          (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                           (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) +
                          (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                           pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2)) /
                              (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          2.5))) /
                             ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                         (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                         (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(y * pi) * cos(y * pi) * sin(x * pi) *
                          sin(x * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                         (4 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                          ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                            (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) -
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)))) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) -
                    (((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                       (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                          (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                     1.5)) -
                      (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                       (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                     1.5))) /
                         ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                     (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                      ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       0.5)) /
                         (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                        ((-2 * x + (14 * sin(t * pi)) / 25 + 1) / (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) +
                                                                   pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                         (pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                             pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                  pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                 2)) -
                    (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                    ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)))) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                    (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                            2)) +
               pi * (x - 0.5) *
                   ((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5))) /
                        ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                    (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                     ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                      (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      0.5)) /
                        (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) +
               pi *
                   (((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                      (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5))) /
                        ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                    (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                     ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                      (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      0.5)) /
                        (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                   (y - 0.5) +
               ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                ((pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) /
                      (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                  1) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3))) -
                 ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                 (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (-2 * x + (14 * sin(t * pi)) / 25 + 1) * pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)), 2) +
                 ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                  (((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                        (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                    (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) +
                    (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                     pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (1000 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             2.5)) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) +
                    (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2)) /
                        (1000 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             2.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                   (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(y * pi) * cos(y * pi) * sin(x * pi) *
                    sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                   (4 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                      (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)))) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)))) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                 ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                  (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         2))) /
                   (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
               ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                ((pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) /
                      (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                  1) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3))) -
                 (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                 ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                  (-((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                         (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) +
                     (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                      pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (1000 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              2.5)) +
                     (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                      pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (1000 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              2.5)) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                   (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(x * pi) * sin(y * pi) *
                    sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                   (4 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)))) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)))) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                 ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (2 * y + (14 * cos(t * pi)) / 25 - 1) * pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)), 2) +
                 (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (2 * y + (14 * cos(t * pi)) / 25 - 1) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                  (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         2))) /
                   (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
               (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                   (250 *
                    sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
               (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                   (250 *
                    sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
               (4 * pi * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                 (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     (250 *
                      sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                 (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                     (250 *
                      sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                 0.5)) /
                   (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2));
    }

    R fun_neumann_left(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
    }

    R fun_neumann_bottom(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
    }

    R fun_neumann_right(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
    }

    R fun_neumann_top(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return -(pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
    }
} 


//* Set numerical example (options: "ex1", "ex2", or "ex3")
#define ex1
//* Set scheme for the method (options: "classical", "conservative")
#define classical
//* Set stabilization method (options: "fullstab", "macro")
#define fullstab

#define use_h    // to set mesh size using the h parameter. Write use_n to decide
                 // using nx, ny.
#define use_tnot // write use_t to control dT manually. Otherwise it is set
                 // proportional to h.
#define conservation

// Setup two-dimensional class types
const int d = 2;

typedef MeshQuad2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

using namespace Example1;


int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   
    int thread_count = 1;
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.1;             // starting mesh size
    double dT = 0.5;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;
    const  double tfinal = .12;

    // Time integration quadrature
    const size_t quadrature_order_time = 9;
    const QuadratureFormular1d &qTime(*Lobatto(quadrature_order_time));  // specify order of quadrature in time
    const Uint nbTime       = qTime.n;
    const Uint lastQuadTime = nbTime - 1;

    // Matrix to hold non-linear part in Newton iteration
    std::vector<std::map<std::pair<int, int>, double>> jacobian(thread_count);

    // Space integration quadrature 
    Levelset<2> phi;
    ProblemOption option;
    const int quadrature_order_space       = 9;
    option.order_space_element_quadrature_ = quadrature_order_space;
    AlgoimCutFEM<Mesh, Levelset<2>> convdiff(qTime, phi, option);

    // Global parameters
    const double tau_F_bulk = 1.;    // face stabilization
    const double tau_F_surf = 1.;    // face stabilization
    const double tau_G = 1.;         // interface stabilization
    
    const double delta_bulk = 0.3;   // macro parameter
    const double delta_surf = 0.5;   // macro parameter

    const double D = 0.01;           // diffusion bulk
    const double DGamma = 1.;        // diffusion surface
    
    std::string ex, method, stab;

    // Paths to store data
    ex = "example1";

    const std::string path_output_data    = "../output_files/langmuir/" + ex + "/data/";
    const std::string path_output_figures = "../output_files/langmuir/" + ex + "/paraview/";

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    // Data file to hold problem data
    #ifdef conservative
    method = "conservative";
    #else
    method = "non_conservative";
    #endif

    #ifdef fullstab
    stab = "full";
    #else
    stab = "macro_dsurf_" + std::to_string(delta_surf);
    #endif

    std::ofstream output_data(path_output_data + "data_" + method + "_" + stab + ".dat", std::ofstream::out);
    
    output_data << method << ",\t";
    output_data << stab << ",\t";
    output_data << "tau_F_bulk = " << tau_F_bulk << ",\t tau_F_surf = " << tau_F_surf << ",\t tau_G = " << tau_G << ",\t N = " << quadrature_order_time << ",\t T = " << tfinal << ",\t Example: " << ex;
    output_data << "\n---------------------\n";
    output_data << "h, \t dt,   L2(Omega(T)), L2(Gamma(T)), L2(Omega(t),0,T), L2(Gamma(t),0,T), e_c(T)\n";
    output_data.flush();

    // Data file to hold DOF indices
    std::ofstream indices(path_output_data + "indices.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors, errors_surf, errors_T, errors_T_surf, hs, nxs, nys, dts, omega, gamma,
        global_conservation_errors, reynold_error;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh

        const double x0 = 0., y0 = 0., lx = 1., ly = 1.;

#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h                          = lx / (nx - 1);
#endif

        const Mesh Th(nx, ny, x0 - Epsilon, y0 - Epsilon, lx, ly);

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 3;

        dT = h / divisionMeshSize;

        total_number_iteration = int(tfinal / dT);
#endif
        dT        = tfinal / total_number_iteration;
        time_step = dT;

        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;
        dts.at(j) = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "Iteration " << j + 1 << "/" << iterations << '\n';
        }

        std::cout << "h  = " << h << '\n';
        std::cout << "nx = " << nx << '\n';
        std::cout << "ny = " << ny << '\n';
        std::cout << "dT = " << dT << '\n';

        FESpace Vh(Th, DataFE<Mesh>::P2); // Background FE Space

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly); // FE Space in time
        const Uint ndf_time_slab      = Ih[0].NbDoF();

        // Velocity field
        LagrangeQuad2 FEvelocity(1);
        FESpace VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Visualization
        FESpace Lh(Th, DataFE<Mesh>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double mass_last_previous, mass_last_previous_surf, mass_initial, mass_initial_surf;
        double intF = 0, int_outflow = 0, intF_surf = 0, intF_total = 0,
               intF_surf_total                = 0; // hold integrals of rhs and Neumann bcs
        double global_conservation_error = 0, local_conservation_error = 0, error_bulk = 0., error_surf = 0.,
               error_I = 0., error_I_surf = 0.;
        std::vector<double> local_conservation_errors;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * time_step;
            const TimeSlab &In(Ih[iter]);
            // const TimeSlab &In_interpolation(Ih_interpolation[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * time_step << '\n';
            std::cout << "dT = " << dT << '\n';

            // initialization of the interface in each quadrature point

            for (int i = 0; i < nbTime; ++i) {

                R tt  = In.Pt(R1(qTime(i).x));
                phi.t = tt;

                interface.init(i, Th, phi);

                ls[i].init(Lh, fun_levelSet, tt);
            }

            // Create active meshes
            ActiveMesh<Mesh> Thi(Th);

            Thi.truncate(interface, 1);

            //  Cut FE space
            CutSpace Wh(Thi, Vh);

            // Surface active mesh
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);

            // Surface FE space
            CutFESpace WhGamma(ThGamma, Vh);

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Add time slab to surface space on top of the data
            convdiff.add(WhGamma, In);

            // Objects needed for the weak form
            Normal n;
            
            // Data for initial solution
            int idx_s0  = Wh.NbDoF() * In.NbDoF();      // index where uS^0 start in the solution array
            int idx_s1  = WhGamma.NbDoF() * In.NbDoF(); // degrees of freedom in surface x time space
            int tot_dof = (int)convdiff.get_nb_dof();   // total degrees of freedom

            std::vector<double> data_init(convdiff.get_nb_dof(), 0.); // initial data total
            std::span<double> data_init_span(data_init);
            std::span<double> data_uhB0 = std::span<double>(data_init.data(), Wh.NbDoF());
            std::span<double> data_uhS0 = std::span<double>(data_init.data() + idx_s0, WhGamma.NbDoF());

            if (iter == 0) {
                interpolate(Wh, data_uhB0, fun_uBulkInit);
                interpolate(WhGamma, data_uhS0, fun_uSurfInit);
            } else {
                convdiff.initialSolution(data_init_span);
            }

            std::vector<double> data_all(data_init);
            std::span<double> data_uhB = std::span<double>(data_all.data(), Wh.NbDoF()*In.NbDoF());
            std::span<double> data_uhS = std::span<double>(data_all.data() + idx_s0, WhGamma.NbDoF()*In.NbDoF());

            Fun_h uhB0(Wh, data_uhB0);
            Fun_h uhS0(WhGamma, data_uhS0);
            Fun_h uhB(Wh, In, data_uhB);
            Fun_h uhS(WhGamma, In, data_uhS);

            int newton_iterations = 0;
            double tol            = 1e-10;
            double residual       = 1.;

            // Test functions
            FunTest uB(Wh, 1), uS(WhGamma, 1), vB(Wh, 1), vS(WhGamma, 1);

            while (1) {

                // Assemble linear part of Jacobian DF
                if (newton_iterations == 0) {
                    
                #ifdef conservative

                     // Time derivative terms
                    convdiff.addBilinear(-innerProduct(uB, dt(vB)), Thi, In);
                    convdiff.addBilinear(-innerProduct(uS, dt(vS)), interface, In);

                    // Time penalty terms (only) LHS
                    convdiff.addBilinear(+innerProduct(uB, vB), Thi, lastQuadTime, In);
                    convdiff.addBilinear(+innerProduct(uS, vS), *interface(lastQuadTime), In, lastQuadTime);

                    // Convection
                    convdiff.addBilinear(-innerProduct(uB, (vel.exprList() * grad(vB))), Thi, In);
                    convdiff.addBilinear(-innerProduct(uS, (vel.exprList() * grad(vS))),
                                         interface, In);


                #elif defined(classical)

                    // Time derivative terms
                    convdiff.addBilinear(+innerProduct(dt(uB), vB), Thi, In);
                    convdiff.addBilinear(+innerProduct(dt(uS), vS), interface, In);

                    // Time penalty terms (only) LHS
                    convdiff.addBilinear(+innerProduct(uB, vB), Thi, 0, In);
                    convdiff.addBilinear(+innerProduct(uS, vS), *interface(0), In, 0);

                    // Convection
                    convdiff.addBilinear(+innerProduct((vel.exprList() * grad(uB)), vB), Thi, In);
                    
                    convdiff.addBilinear(
                        + innerProduct((vel.exprList() * grad(uS)), vS) 
                        + innerProduct(uS * divS(vel), vS)
                        , interface
                        , In);
                #endif

                    // Diffusion
                    convdiff.addBilinear(+innerProduct(D * grad(uB), grad(vB)), Thi, In);
                    convdiff.addBilinear(+innerProduct(DGamma * gradS(uS), gradS(vS)), interface, In);

                    // Linear part of coupling term
                    convdiff.addBilinear(innerProduct(uB - uS, vB - vS), interface, In);


                #if defined(fullstab)

                    convdiff.addPatchStabilization(
                        + innerProduct(tau_F_bulk / h / h * jump(uB), jump(vB))
                        , Thi
                        , In
                    );

                    convdiff.addPatchStabilization(
                        + innerProduct(tau_F_surf / h / h / h * jump(uS), jump(vS))
                        , ThGamma
                        , In
                    );

                    // convdiff.addFaceStabilization(
                    //     + innerProduct(tau_F_bulk * h * jump(grad(uB) * n), jump(grad(vB) * n)) 
                    //     + innerProduct(tau_F_bulk * h * h * h * jump(grad(grad(uB) * n) * n), jump(grad(grad(vB) * n) * n)),
                    //     Thi, In);

                    // convdiff.addFaceStabilization(
                    //     + innerProduct(tau_F_surf * jump(grad(uS) * n), jump(grad(vS) * n)) 
                    //     + innerProduct(tau_F_surf * h * h * jump(grad(grad(uS) * n) * n), jump(grad(grad(vS) * n) * n)),
                    //     ThGamma, In);


                #elif defined(macro)
                
                    AlgoimMacro<Mesh, Levelset<2>> TimeMacro(Thi, 0.4, phi, In, qTime);
                    TimeMacro.findSmallElement();
                    TimeMacro.createMacroElement();
                    TimeMacro.setInnerEdges();

                    AlgoimMacroSurface<Mesh, Levelset<2>> TimeMacroSurf(ThGamma, 0.5, phi, In, qTime);
                    TimeMacroSurf.findSmallElement();
                    TimeMacroSurf.createMacroElement();
                    TimeMacroSurf.setInnerEdges();

                    convdiff.addPatchStabilization(
                        + innerProduct(tau_F_bulk / h / h * jump(u), jump(v))
                        , Thi
                        , In
                        , TimeMacro
                    );

                    convdiff.addPatchStabilization(
                        + innerProduct(tau_F_surf / h / h / h * jump(uS), jump(vS))
                        , ThGamma
                        , In
                        , TimeMacroSurf
                    );

                #endif

                    convdiff.addBilinear(
                        + innerProduct(tau_G * grad(uS) * n, grad(vS) * n)
                        + innerProduct(tau_G * h * h * grad(grad(uS) * n) * n, grad(grad(vS) * n) * n)
                        , interface
                        , In);

                }

                //convdiff.gather_map();  //? should this be here??

                // Multiply the linear part of the Jacobian by the initial guess to evaluate in the initial guess (which
                // is updated in each Newton iteration)
                convdiff.addMatMul(data_all);   // -> goes into the rhs (this does not change the matrix convdiff.mat_)

                // Add the remaining terms of F (the residual)
                convdiff.addLinear(- innerProduct(uhB.expr() * uhS.expr(), vB - vS), interface, In);

                // Time penalty RHS
                if (iter == 0) {    
                    convdiff.addLinearExact(fun_uBulk, -innerProduct(1., vB), Thi, 0, In);
                    convdiff.addLinearExact(fun_uSurf, -innerProduct(1., vS), *interface(0), In, 0);

                }
                else {
                    convdiff.addLinear(-innerProduct(uhB0.expr(), vB), Thi, 0, In);
                    convdiff.addLinear(-innerProduct(uhS0.expr(), vS), *interface(0), In, 0);
                }

                // RHS force functions
                convdiff.addLinearExact(fun_rhsBulk, -innerProduct(1., vB), Thi, In);
                convdiff.addLinearExact(fun_rhsSurf, -innerProduct(1., vS), interface, In);

                // Switch the system matrix to the jacobian instead of convdiff.mat_
                convdiff.set_map(jacobian);

                // Initialize the Jacobian as the convdiff matrix since these terms are also in the Jacobian
                jacobian[0] = convdiff.mat_[0];

                // Assemble remaining non-linear part of Jacobian (evaluated in the initial guess so it becomes linear)
                convdiff.addBilinear(
                    - innerProduct(uB * uhS.expr(), vB - vS)
                    - innerProduct(uS * uhB.expr(), vB - vS)
                    , interface
                    , In);
                           

                // Solve DF(u0)(w) = F(u0)
                convdiff.solve(jacobian[0], convdiff.rhs_);

                // Retrieve the bulk solution
                std::span<double> dwB = std::span<double>(convdiff.rhs_.data(), Wh.NbDoF() * In.NbDoF());
                std::span<double> dwS = std::span<double>(convdiff.rhs_.data() + idx_s0, WhGamma.NbDoF() * In.NbDoF());

                // Compute the residuals in max norms
                // const auto max_res_B = std::max_element(dwB.begin(), dwB.end());
                // const auto max_res_S = std::max_element(dwS.begin(), dwS.end());

                // Compute the residual as the maximum of the absolute values of the residual data
                const auto max_res_B = std::max_element(dwB.begin(), dwB.end(), [](int a, int b) { return std::abs(a) < std::abs(b); });
                const auto max_res_S = std::max_element(dwS.begin(), dwS.end(), [](int a, int b) { return std::abs(a) < std::abs(b); });
                
                residual = std::max(*max_res_B, *max_res_S);

                std::cout << " Residual \t : \t" << residual << '\n';

                // Update the solution
                std::span<double> dw = std::span<double>(convdiff.rhs_.data(), convdiff.get_nb_dof());
                std::transform(data_all.begin(), data_all.end(), dw.begin(), data_all.begin(),
                               [](double a, double b) { return a - b; });

                // Reset the rhs vector
                convdiff.rhs_.resize(convdiff.get_nb_dof());
                convdiff.rhs_ = 0.0; 

                newton_iterations++;

                if (std::abs(residual) < tol) {
                    std::span<double> data_all_span(data_all);
                    convdiff.saveSolution(data_all_span);
                    convdiff.cleanBuildInMatrix();
                    convdiff.set_map();

                    break;
                }

                if (newton_iterations >= 4) {
                    std::cout << "Newton's method didn't converge in " + std::to_string(newton_iterations) + " iterations, breaking loop. \n";
                    break;
                }

            }

            // Compute L2(Omega(t), 0, T) and L2(Gamma(t), 0, T)     
            Fun_h fun_uhB_t(Wh, In, data_uhB);
            Fun_h fun_uhS_t(WhGamma, In, data_uhS);
            error_I += L2_norm_T(fun_uhB_t, fun_uBulkD, Thi, In, qTime, phi, quadrature_order_space);
            error_I_surf += L2_norm_surf_T(fun_uhS_t, fun_uSurf, interface, In, qTime, phi, quadrature_order_space); 

            // Compute area of domain in time quadrature point 0
            Fun_h funone(Wh, fun_one);

            double intGamma = integral_algoim(funone, In, interface, phi, 0, quadrature_order_space) / dT;
            double intOmega = integral_algoim(funone, Thi, phi, In, qTime, lastQuadTime, quadrature_order_space);

            gamma[j] = intGamma;
            omega[j] = intOmega;


            std::vector<double> sol_uhB(Wh.get_nb_dof());
            std::vector<double> sol_uhS(WhGamma.get_nb_dof()); 

            for (int n = 0; n < ndf_time_slab; ++n) {
                // get the DOFs of u corresponding to DOF n in time and sum with the
                // previous n
                std::vector<double> uB_dof_n(data_uhB.begin() + n * Wh.get_nb_dof(),
                                            data_uhB.begin() + (n + 1) * Wh.get_nb_dof());
                std::transform(sol_uhB.begin(), sol_uhB.end(), uB_dof_n.begin(), sol_uhB.begin(), std::plus<double>());

                // get the DOFs of p corresponding to DOF n in time and sum with the
                // previous n
                std::vector<double> uS_dof_n(data_uhS.begin() + n * WhGamma.get_nb_dof(),
                                            data_uhS.begin() + (n + 1) * WhGamma.get_nb_dof());
                std::transform(sol_uhS.begin(), sol_uhS.end(), uS_dof_n.begin(), sol_uhS.begin(), std::plus<double>());
            }            

            Fun_h fun_uhB(Wh, sol_uhB);
            Fun_h fun_uhS(WhGamma, sol_uhS);

            error_bulk = L2_norm_cut(fun_uhB, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1, quadrature_order_space);
            error_surf = L2_norm_surface(fun_uhS, fun_uSurf, *interface(lastQuadTime), current_time + dT, phi, 0, 1);

            std::cout << " t_n -> || u-uex||_Omega(T) = " << error_bulk << '\n';
            std::cout << " t_n -> || u-uex||_Gamma(T) = " << error_surf << '\n';

            errors[j]      = error_bulk;
            errors_surf[j] = error_surf;

// Compute conservation error
// if (iterations == 1 && h > 0.01) {
#if defined(conservation)
            intF = integral_algoim(fun_rhsBulk, 0, Thi, phi, In, qTime,
                                   quadrature_order_space); // integrate source over In
            intF_surf = integral_algoim(fun_rhsSurf, In, interface, phi, 0,
                                   quadrature_order_space); // integrate flux boundary over In

            intF_total += intF;
            intF_surf_total += intF_surf;

            double mass_last = integral_algoim(fun_uhB, Thi, phi, In, qTime, lastQuadTime,
                                               quadrature_order_space); // mass in last quad point

            double mass_last_surf = integral_algoim(fun_uhS, *interface(lastQuadTime), 0, phi, In, qTime,
                                                    lastQuadTime, quadrature_order_space);

            if (iter == 0) {
                mass_initial       = integral_algoim(fun_uBulk, Thi, phi, In, qTime, 0, quadrature_order_space);
                mass_initial_surf = integral_algoim(fun_uSurf, *interface(0), 0, phi, In, qTime,
                                                    0, quadrature_order_space);
                mass_last_previous = mass_initial;
                mass_last_previous_surf = mass_initial_surf;
                // mass_last_previous = integral_algoim(b0h, Thi, phi, In, qTime, 0);
            }

            local_conservation_error  = (mass_last + mass_last_surf - mass_last_previous - mass_last_previous_surf - intF - intF_surf);
            // global_conservation_error += local_conservation_error;
            global_conservation_error = (mass_last + mass_last_surf - mass_initial - mass_initial_surf - intF_total - intF_surf_total);

            std::cout << "global_conservation_error: " << global_conservation_error << "\n";
            std::cout << "local_conservation_error: " << local_conservation_error << "\n";

            // output_data << std::setprecision(10);
            // output_data << current_time << "," << (mass_last - mass_last_previous) << "," << intF << "," << intG <<
            // ","
            //            << local_conservation_error << '\n';

            mass_last_previous = mass_last; // set current last to previous last for next time slab
                                            //}
            mass_last_previous_surf = mass_last_surf;

            global_conservation_errors[j] = std::fabs(global_conservation_error);
            local_conservation_errors.push_back(std::fabs(local_conservation_error));

#endif

            if ((iterations == 1) && (h > 0.01)) {
                
                Paraview<Mesh> writerTh(Th, path_output_figures + "Th.vtk");
                writerTh.add(ls[0], "levelSet0", 0, 1);
                writerTh.add(ls[1], "levelSet1", 0, 1);
                writerTh.add(ls[2], "levelSet2", 0, 1);
                
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_" + std::to_string(iter + 1) + ".vtk");
                Paraview<Mesh> writer_surf(ThGamma, path_output_figures + "surf_" + std::to_string(iter + 1) + ".vtk");
                // writer.add(b0h, "bulk", 0, 1);
                // writer.add(sol_h, "bulk_end", 0, 1);

                writer.add(fun_uhB, "bulk", 0, 1);
                writer_surf.add(fun_uhS, "surf", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                
                Fun_h uSex(WhGamma, fun_uSurf, current_time);
                Fun_h fS(WhGamma, fun_rhsSurf, current_time);

                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);

                writer_surf.add(uSex, "surf_exact", 0, 1);
                writer_surf.add(fS, "surf_rhs", 0, 1);

                //writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
                writer.add(fabs(fun_uhB.expr() - uBex.expr()), "bulk_error");
                writer_surf.add(fabs(fun_uhS.expr() - uSex.expr()), "surf_error");
                
                
                writer.writeActiveMesh(Thi, path_output_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writer.writeFaceStab(Thi, 0, path_output_figures + "Edges" + std::to_string(iter + 1) + ".vtk");

                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 0, 2,
                                             path_output_figures + "AlgoimQuadrature_0_" + std::to_string(iter + 1) +
                                                 ".vtk");
                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 1, 2,
                                             path_output_figures + "AlgoimQuadrature_1_" + std::to_string(iter + 1) +
                                                 ".vtk");
                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 2, 2,
                                             path_output_figures + "AlgoimQuadrature_2_" + std::to_string(iter + 1) +
                                                 ".vtk");
            }

            // if (iterations > 1 && iter == total_number_iteration - 1)
            //     output_data << h << "," << dT << "," << errBulk << '\n';

            iter++;
        }

        errors_T[j] = std::sqrt(error_I);
        errors_T_surf[j] = std::sqrt(error_I_surf);

        output_data << h << "," << dT << "," << error_bulk << "," << error_surf << "," << errors_T[j] << "," << errors_T_surf[j] << "\n";
        output_data.flush();

        std::cout << "error_T = " << errors_T[j] << "\n";
        std::cout << "error_T_surf = " << errors_T_surf[j] << "\n";

        // std::cout << "\n";
        // std::cout << "Local conservation error = [";
        // for (auto &err : local_conservation_errors) {

        //     std::cout << err;

        //     std::cout << ", ";
        // }
        // std::cout << "]\n";

// Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
#elif defined(use_t)
        dT *= 0.5;
#elif defined(use_h)
        h *= 0.5;
        // h *= sqrt(0.5);
#endif
    }

    std::cout << std::setprecision(16);
    std::cout << '\n';
    std::cout << "Errors Bulk = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors Surface = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_surf.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors In= [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_T.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors In Surface = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_T_surf.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Global Conservation Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << global_conservation_errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Reynold error = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << reynold_error.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "|Gamma_h| = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << gamma.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "|Omega_h| = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << omega.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "dT = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << dts.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << '\n';

    std::cout << "nx = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << nxs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << '\n';

    std::cout << "ny = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << nys.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    return 0;
}
