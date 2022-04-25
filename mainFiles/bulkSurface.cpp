#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../util/redirectOutput.hpp"
#include "../num/matlab.hpp"
#include "paraview.hpp"


namespace TestSaraExample3 {
  R fun_levelSet(const R2 P, const int i) { return sqrt((P.x-0.1)*(P.x-0.1) + (P.y-0.0)*(P.y-0.0)) - 0.3;}
  R fun_rhs(const R2 P,const int i) { R2 R(0,-0.98); return (i<2)?R[i] : 0;}
  R fun_velocity(const R2 P, const int i){
    if(i == 0) return -0.5*(1+cos(M_PI*P.x))*sin(M_PI*P.y);
    else       return  0.5*(1+cos(M_PI*P.y))*sin(M_PI*P.x);
  }
  R fun_init_surface(const R2 P, int i) { return 0.;}
  R fun_w(double r) {
    return 0.5*(1-cos((r-0.3)*M_PI/(0.5*0.3)));
  }
  R fun_init_bulk(const R2 P, int i) {
    double r = sqrt((P.x-0.1)*(P.x-0.1) + P.y*P.y);
    if(r > 1.5*0.3) return 0.5*(1-P.x*P.x)*(1-P.x*P.x);
    else if(0.3 < r && r <= 1.5*0.3) return 0.5*(1-P.x*P.x)*(1-P.x*P.x)*fun_w(r);
    else return 0.;}

  R sigma(const R w) {return 24.5 - 0.5*w;}
  R Dsigma(const R w) {return - 0.5*w;}

  R fun_x(const R2 P,const int i) {return P[i];}
  R fun_1(const R2 P,const int i) {return 1;}
}

namespace TestSaraExample1 {

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
        return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.22)*(P.y-0.22)) - 0.17;
    }

    // double fun_levelSet(const R2 P, const int i, const R t) {
    //     R xc = 0.5, yc = 0.5;
    //     return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
    // }
    //
    // // Level-set function initial
    // double fun_levelSet(const R2 P, const int i) {
    //     return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.17;
    // }

    // The rhs Neumann boundary condition
    R fun_rhs_neumann_bc(const R2 P, const int i) {return 0.;}

    // Velocity field
    R fun_velocity(const R2 P, const int i) {
        if (i == 0) return M_PI*(0.5-P.y);
        else return M_PI*(P.x-0.5);
    }

    // Normal x-direction
    R n1(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return -(P.x - xc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Normal y-direction
    R n2(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return -(P.y - yc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Initial solution bulk
    R fun_init_bulk(const R2 P, const int i) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y);
        return P.y/2+2;

    }

    // Initial solution surface
    R fun_init_surface(const R2 P, const int i) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)
        // + M_PI/250*sin(M_PI*P.x)*cos(M_PI*P.y)*n1(P, 0)
        // + M_PI/250*cos(M_PI*P.x)*sin(M_PI*P.y)*n2(P, 0);
        return P.y/2+2;
    }

    // Exact solution bulk
    R fun_bulk(const R2 P, const int i, const R t) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
        return P.y/2+2;

    }
    R fun_gradbulk(const R2 P, const int i, const R t) {
        if(i==0) return -0.4*M_PI*sin(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
        else return -0.4*M_PI*cos(M_PI*P.x)*sin(M_PI*P.y)*cos(2*M_PI*t);
    }

    R fun_bulk2(const R2 P, const int i, const int d, const R t) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
        return P.y/2+2;

    }

    // Exact solution surface
    R fun_surface(const R2 P, const int i, const R t) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t)
        // + M_PI/250*sin(M_PI*P.x)*cos(M_PI*P.y)*n1(P, t)
        // + M_PI/250*cos(M_PI*P.x)*sin(M_PI*P.y)*n2(P, t);
        return P.y/2+2;
    }

    R fun_rhs_bulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;
        return   (pi*(x - 1/2))/2;
    }



    R fun_rhs_surf(const R2 P, const int i, const R t) {
      R x = P.x, y = P.y;
      return (1250*y + 350*cos(M_PI*t) - 625)/(700*sin(M_PI*t) - 2500*y - 700*cos(M_PI*t) - 2500*x + 1400*y*cos(M_PI*t) - 1400*x*sin(M_PI*t) + 2500*x*x + 2500*y*y + 1446) + (M_PI*(x - 1./2))/2;
    }

    // RHS fB bulk
    // R fun_rhs_bulk(const R2 P, const int i, const R t) {
    //     R x = P.x, y = P.y;
    //     return (M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*cos(M_PI*y))/125
    //     - (4*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(2*M_PI*t))/5
    //     - (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*sin(M_PI*y)*(x - 0.5))/5
    //     + (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*y)*sin(M_PI*x)*(y - 0.5))/5;
    //
    // }

    // RHS fS surface
    // R fun_rhs_surf(const R2 P, const int i, const R t) {
    //     R x = P.x, y = P.y;
    //
    //     return -(M_PI*(46875*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 46875*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 2450*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 7350*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 1372*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 31250*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 15625*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 15625*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 93750*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 31250*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 31250*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 3125000*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 3125000*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 31250*x*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 125000*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 125000*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 31250*y*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 35000*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 8750*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 8750*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 35000*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 31250*x*x*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 93750*x*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*x*x*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 93750*y*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 31250*y*y*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*y*y*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 7350*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 2450*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 1372*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 1750000*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 3500000*cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 52500*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 17500*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 35000*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 17500*x*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 17500*x*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 52500*y*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 17500*y*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 1750000*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3500000*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 17500*x*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 52500*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 17500*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 17500*y*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 62500*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 62500*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 4900*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 4900*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 4900*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 4375*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 31250*x*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 31250*x*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 93750*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 31250*x*x*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*x*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 31250*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 31250*y*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 31250*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 93750*y*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 62500*y*y*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 4900*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) + 2450*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 7350*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 8750*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 8750*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 4375*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 1372*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 31250*x*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 62500*x*x*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 187500*x*x*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 125000*x*x*x*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 31250*y*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 187500*y*y*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 62500*y*y*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 125000*y*y*y*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 12500000*x*x*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 12500000*y*y*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 7350*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 2450*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 14700*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 4900*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 4375*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 8750*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 8750*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 1372*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 2744*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 62500*x*y*y*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*x*x*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 980000*cos(t*M_PI)*cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 4900*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 14700*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 4375*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 2744*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 4900*x*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 17500*x*x*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 14700*y*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 52500*y*y*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 980000*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 14700*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 52500*x*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 4900*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 17500*y*y*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 1372*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 1372*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 2450*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 686*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 31250*x*x*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 93750*x*x*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*x*x*x*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 93750*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 31250*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*y*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 2450*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) + 7350*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 2450*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 686*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI) - 1372*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 93750*x*x*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 125000*x*x*x*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 62500*x*x*x*x*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 93750*y*y*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 125000*y*y*y*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 62500*y*y*y*y*M_PI*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 2450*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 7350*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 2450*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 1372*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 686*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 3125000*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 2450*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 686*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 62500*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 62500*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 1562500*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 1562500*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 17500*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 62500*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 31250*x*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 31250*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 62500*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 3125000*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 8750*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 17500*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 17500*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 187500*x*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 62500*x*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 125000*x*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 187500*y*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 125000*y*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 62500*y*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 6250000*x*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 12500000*x*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 12500000*y*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 52500*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 35000*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 17500*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 17500*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 8750*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 62500*x*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 62500*x*y*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 4900*x*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 26250*x*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 1372*x*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 17500*x*x*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 4900*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 8750*y*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 62500*x*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 62500*x*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 4900*x*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) + 4900*x*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 17500*x*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 8750*x*x*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 4900*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) - 14700*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 52500*y*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 26250*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 1372*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI) - 17500*y*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 686*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) + 686*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 14700*x*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 4900*x*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 8750*x*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 52500*x*x*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 4900*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 14700*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 78750*y*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 17500*y*y*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 1372*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 52500*y*y*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 6250000*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 1372*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 1372*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 14700*x*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 78750*x*x*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 1372*x*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 52500*x*x*x*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 4900*y*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 8750*y*y*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 1750000*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 3125000*x*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3125000*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 6250000*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 686*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 686*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 35000*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 62500*x*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 62500*x*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 875000*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 875000*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 1750000*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*x*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*y*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 17500*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 17500*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 35000*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 17500*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 52500*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 375000*x*y*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 125000*x*y*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 125000*x*y*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 1750000*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 875000*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 875000*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 105000*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 35000*x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 35000*x*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 52500*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 17500*x*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 105000*y*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 35000*y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 17500*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 17500*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 1750000*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 7000000*y*cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 4900*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 4900*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 35000*x*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 105000*x*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 105000*y*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 35000*y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 35000*y*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 35000*x*y*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 7000000*x*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 4900*x*x*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 14700*y*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 29400*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 9800*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 9800*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 6250000*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 6250000*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 35000*x*y*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 14700*x*x*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 4900*y*y*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 490000*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 9375000*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3125000*x*x*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*x*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 3125000*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 9375000*y*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 6250000*y*y*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 9800*x*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 9800*y*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 17500*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 8750*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 62500*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 62500*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 62500*x*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) + 62500*x*x*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 490000*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 245000*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 245000*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 4900*x*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 17500*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 17500*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 8750*x*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 17500*x*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 14700*y*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 52500*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 17500*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) - 17500*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 52500*y*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 125000*x*y*y*M_PI*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 125000*x*x*y*M_PI*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 245000*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 245000*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 14700*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 9800*x*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 8750*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 17500*x*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 52500*x*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 35000*x*x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 52500*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 4900*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 29400*y*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 35000*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 17500*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 17500*y*M_PI*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 105000*y*y*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 17500*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 1372*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 1372*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 4900*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) + 4900*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 29400*x*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 35000*x*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 105000*x*x*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 9800*y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 8750*y*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 35000*y*y*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 2744*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 2744*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 19600*x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) - 19600*y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) + 6250000*x*y*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 6250000*x*x*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 17500*x*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 490000*x*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 490000*y*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 3500000*y*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 35000*x*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI) + 17500*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 490000*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 3500000*x*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 490000*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 9800*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 9800*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 17500*x*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 35000*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 9800*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 9800*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) - 17500*x*y*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 9800*x*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) - 17500*x*y*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI) + 9800*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) - 17500*x*x*y*M_PI*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 1372*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI) + 9800*x*x*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 1372*y*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) - 9800*y*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI) + 17500*x*x*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI) + 17500*x*y*y*M_PI*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 3500000*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 6250000*x*y*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 6250000*x*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 1372*x*M_PI*M_PI*cos(t*M_PI)*cos(t*M_PI)*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 1372*y*M_PI*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI) - 1750000*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3500000*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 1750000*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3500000*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 12500000*x*y*M_PI*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 35000*x*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI) - 3500000*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3500000*x*M_PI*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 1750000*x*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 1750000*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 70000*x*y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(2*t*M_PI)*sin(x*M_PI) - 35000*x*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 3500000*y*M_PI*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) - 9800*x*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI) + 9800*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI) - 70000*x*y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(2*t*M_PI)*sin(y*M_PI) + 980000*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3500000*x*y*M_PI*cos(t*M_PI)*cos(2*t*M_PI)*cos(x*M_PI)*sin(y*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5)))) + 3500000*x*y*M_PI*cos(2*t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt((((y + (7*cos(t*M_PI))/25 - 0.5)*(y + (7*cos(t*M_PI))/25 - 0.5)) + ((- x + (7*sin(t*M_PI))/25 + 0.5)*(- x + (7*sin(t*M_PI))/25 + 0.5))))))/(12500*sqrt(((((7*sin(M_PI*t))/25 - x + 0.5)*((7*sin(M_PI*t))/25 - x + 0.5)) + ((y + (7*cos(M_PI*t))/25 - 0.5)*(y + (7*cos(M_PI*t))/25 - 0.5))))*(- 1250*x - 1250*y - 350*cos(t*M_PI) + 350*sin(t*M_PI) + 700*y*cos(t*M_PI) - 700*x*sin(t*M_PI) + 1250*x*x + 1250*y*y + 98*cos(t*M_PI)*cos(t*M_PI) + 98*sin(t*M_PI)*sin(t*M_PI) + 625));
    // }

    // Set rhs Neumann boundary condition (they are zero on the boundary)
    R fun_g_bottom(const R2 P, const int i, const R t) {
        return -(M_PI*cos(2*M_PI*t)*cos(M_PI*P.x)*sin(M_PI*P.y))/250;
    }

    R fun_g_right(const R2 P, const int i, const R t) {
        return (M_PI*cos(2*M_PI*t)*cos(M_PI*P.y)*sin(M_PI*P.x))/250;
    }

    R fun_g_top(const R2 P, const int i, const R t) {
        return (M_PI*cos(2*M_PI*t)*cos(M_PI*P.x)*sin(M_PI*P.y))/250;
    }
    R fun_g_left(const R2 P, const int i, const R t) {
        return -(M_PI*cos(2*M_PI*t)*cos(M_PI*P.y)*sin(M_PI*P.x))/250;
    }
}
namespace TestSaraExample2 {

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        R xc = 0.5 , yc = 0.5;
        return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
      R xc = 0.5 , yc = 0.5;
      return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
     }


    // The rhs Neumann boundary condition
    R fun_rhs_neumann_bc(const R2 P, const int i) {return 0.;}

    // Velocity field
    R fun_velocity(const R2 P, const int i) {
        if (i == 0) return M_PI*(0.5-P.y);
        else return M_PI*(P.x-0.5);
    }

    // Normal x-direction
    R n1(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return -(P.x - xc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Normal y-direction
    R n2(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return -(P.y - yc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Initial solution bulk
    R fun_init_bulk(const R2 P, const int i) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y);
        return P.y/2+2;

    }

    // Initial solution surface
    R fun_init_surface(const R2 P, const int i) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)
        // + M_PI/250*sin(M_PI*P.x)*cos(M_PI*P.y)*n1(P, 0)
        // + M_PI/250*cos(M_PI*P.x)*sin(M_PI*P.y)*n2(P, 0);
        return P.y/2+2;
    }

    // Exact solution bulk
    R fun_bulk(const R2 P, const int i, const R t) {
        // return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
        return P.y/2+2;

    }
    R fun_gradbulk(const R2 P, const int i, const R t) {
        if(i==0) return 0;
        else     return 0.5;
    }

    R fun_bulk2(const R2 P, const int i, const int d, const R t) {
        return P.y/2+2;

    }

    // Exact solution surface
    R fun_surface(const R2 P, const int i, const R t) {
        return P.y/2+2;
    }

    R fun_rhs_bulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;
        return   (pi*(x - 1./2))/2;
    }



    R fun_rhs_surf(const R2 P, const int i, const R t) {
      R x = P.x, y = P.y;
      return (2*y - 1)/(2*(2*x*x - 2*x + 2*y*y - 2*y + 1)) + (pi*(x - 1./2))/2;
    }



    // Set rhs Neumann boundary condition (they are zero on the boundary)
    R fun_g_bottom(const R2 P, const int i, const R t) {
        return -(M_PI*cos(2*M_PI*t)*cos(M_PI*P.x)*sin(M_PI*P.y))/250;
    }

    R fun_g_right(const R2 P, const int i, const R t) {
        return (M_PI*cos(2*M_PI*t)*cos(M_PI*P.y)*sin(M_PI*P.x))/250;
    }

    R fun_g_top(const R2 P, const int i, const R t) {
        return (M_PI*cos(2*M_PI*t)*cos(M_PI*P.x)*sin(M_PI*P.y))/250;
    }
    R fun_g_left(const R2 P, const int i, const R t) {
        return -(M_PI*cos(2*M_PI*t)*cos(M_PI*P.y)*sin(M_PI*P.x))/250;
    }
}


using namespace TestSaraExample2;
int main(int argc, char** argv ){

  const int d = 2;
  typedef Mesh2 Mesh;
  typedef FESpace2 Space;
  typedef CutFESpace<Mesh> CutSpace;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;



  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 20;
  int ny = 20;

  // Mesh2 Kh(nx, ny, -1., -1., 2., 2.);
  Mesh2 Kh(nx, ny, 0., 0., 1., 1.);

  Space Vh(Kh, DataFE<Mesh2>::P1);

  double meshSize(1./(nx-1));
  double hi = 1./(nx-1);
  // Mesh2 Kh("../mesh/RBadaptMesh1_3060.msh");
  // int nx = 30;
  // int ny = 60;
  // double meshSize((1*2)/sqrt(nx*ny));

  int divisionMeshsize = 2;
  double Tend = 0.5;
  double dT = meshSize/divisionMeshsize;
  GTime::total_number_iteration = Tend/dT;
  dT = Tend / GTime::total_number_iteration;
  GTime::time_step = dT;

  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  // Create files to save results
  // std::string meshstr = to_string(nx)+to_string(ny);
  // std::string pathOutpuFolder = "";//../../outputFiles/bulkSurface/example2Sara/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  // std::string pathOutpuFigure = "";//"../../outputFiles/bulkSurface/example2Sara/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  // CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  // myCout << " path to the output files : "+pathOutpuFolder << std::endl;
  // myCout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  // if(MPIcf::IamMaster()) {
  //   std::experimental::filesystem::create_directories(pathOutpuFigure);
  // }
  // std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  // myCout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  // myCout << " Creating the file \"output.txt\" to save cout" << std::endl;
  //
  // myCout << " Mesh size \t" << meshSize << std::endl;
  // myCout << " Number of node \t" << Kh.nv << std::endl;
  // myCout << " Number of element \t" << Kh.nt << std::endl;
  //
  // myCout << " Time Step \t" << dT << std::endl;
  // myCout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  // myCout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  // myCout << " number of quadrature points in time : \t" << nbTime << std::endl;
  // myCout << " number of dof in time per time slab : \t" << ndfTime << std::endl;



  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveSurfactantVTK = true;
  const int frequencyPlottingTheSolution = 2;


  // // Space for the velocity field
  // // ----------------------------------------------
  // myCout << "\n --------------------------------------- \n " << std::endl;
  // myCout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange2 FEvelocity(2);
  Space VelVh(Kh, FEvelocity);
  Fun_h vel(VelVh, fun_velocity);

  // levelSet stuff
  // // ----------------------------------------------
  // myCout << "\n --------------------------------------- \n " << std::endl;
  // myCout << " We use a P2 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;
  Space Lh  (Kh, DataFE<Mesh2>::P1);
  // FESpace2 Lh_k(Kh, DataFE<Mesh2>::P2);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls_k(nbTime), ls(nbTime);
  // for(int i=0;i<nbTime;++i) ls_k[i].init(Lh_k, fun_levelSet);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet);
  // projection(ls_k[0], ls[nbTime-1]);

  // CReinitialization<Mesh2> reinitialization;
  // reinitialization.number_iteration = 2;
  // reinitialization.epsilon_diffusion = 1e-3;
  // reinitialization.dt = dT/4;
  // reinitialization.ON_OFF = "OFF";
  // reinitialization.ODE_method = "Euler";
  // reinitialization.mass_correction = "OFF";
  // reinitialization.precision_correction = 1e-6;
  // reinitialization.max_iteration_correction = 10;
  // reinitialization.info();
  // const int frequencyReinitialization = 5;


  // Declaration of the interface
  TimeInterface<Mesh> interface(qTime);

  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh> bulkSurface(qTime);

  const R epsilon_surfactant = 1.;
  const R Bb = 1.;
  const R Bs = 1.;
  const R Bbs = 1.;
  const R Kb = 0.01;
  const R Ks = 1.;
  const R Kbs = 0.;


  // const Parameter& lambdaG(Parameter::lambdaG);
  // const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  // const CutFEM_Parameter& h(Parameter::h);
  // const CutFEM_Parameter& h_E(Parameter::meas);

  // myCout << " \n Beginning of the time iteration \n"
  //        << " --------------------------------------- \n " << std::endl;


  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {
    double tbegin_iter = CPUtime();
    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter+1 << " / " << GTime::total_number_iteration << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;


    double t_ls0 = CPUtime();
    // ls.begin()->swap(ls[lastQuadTime]);
    for(int i=0;i<nbTime;++i) {

      // projection(ls_k[i], ls[i]);
      // if(iter%frequencyReinitialization == 0 && i == 1) {
      //   reinitialization.perform(ls_k[i], ls[i]);
      // }
      R tt = In.Pt(R1(qTime(i).x));
      ls[i].init(Lh, fun_levelSet, tt);
      interface.init(i,Kh,ls[i]);



      // if(i<lastQuadTime) {
      //   LevelSet::move(ls[i], vel, vel, dt_levelSet, ls[i+1]);
      // }
    }


    ActiveMesh<Mesh> Kh1(Kh);
    Kh1.truncate(interface, -1);
    ActiveMesh<Mesh> Kh0(Kh);
    Kh0.createSurfaceMesh(interface);
    CutSpace Wh(Kh1, Vh);
    CutSpace Sh(Kh0, Vh);

    /*
    PROBLEM DEFINITION
    */

    Fun_h gradUb(VelVh, In, fun_gradbulk);
    Fun_h Us_ex(Sh, In, fun_surface);


    Fun_h fBh(Wh, In, fun_rhs_bulk);          // create a FE-function for the RHS function in the bulk
    Fun_h fSh(Sh, In, fun_rhs_surf);

    Normal n;
    FunTest ub(Wh,1), vb(Wh,1);           // bulk function
    FunTest us(Sh,1), vs(Sh,1);     // surface function

    // Connect spaces to problem
    bulkSurface.initSpace(Wh, In);
    bulkSurface.add(Sh, In);


    // initialize the solution
    int idx_s0 = Wh.NbDoF()*In.NbDoF();
    Rn data_u0(bulkSurface.get_nb_dof(), 0.);
    Rn_ data_b0(data_u0(SubArray(Wh.NbDoF()   ,0)));
    Rn_ data_s0(data_u0(SubArray(Sh.NbDoF(),idx_s0)));

    if(iter == 0) {
      interpolate(Wh, data_b0, fun_init_bulk);
      interpolate(Sh, data_s0, fun_init_surface);
    } else {
      bulkSurface.initialSolution(data_u0);
    }
    Fun_h b0h(Wh, data_b0);
    Fun_h s0h(Sh, data_s0);

    // set uh for newton
    Rn data_uh(data_u0);
    Rn_ data_bh(data_uh(SubArray(Wh.NbDoF()*In.NbDoF(),0)));
    Rn_ data_sh(data_uh(SubArray(Sh.NbDoF()*In.NbDoF(),idx_s0)));
    Fun_h bh(Wh, In, data_bh);
    Fun_h sh(Sh, In, data_sh);
    //
    //
    // int iterNewton  = 0;
    // while(1) {
    //   /*
    //   CONSTRUCTION BILINEAR PART
    //   */
    //   // std::cout << "\n Begin assembly " << std::endl;
    //   double tt0 = CPUtime();
    //   if(iterNewton == 0){

    bulkSurface.addBilinear(
      innerProduct(Bb*dt(ub), vb)
      +innerProduct((vel.expression()*grad(ub)), vb*Bb)
      +innerProduct(grad(ub), grad(vb))*Bb*Kb
      , Kh1
      , In
    );


    bulkSurface.addBilinear(
      innerProduct( ub , Kb*Bb*vb)
      , interface
      , In
    );
    bulkSurface.addLinear(
      innerProduct( Us_ex.expression() , Kb*Bs*vb)
      , interface
      , In
    );
bulkSurface.addLinear(
      innerProduct( gradUb.expression() , Kb*Bs*vb*n)
    , Kh1
    , boundary
    , In
);

bulkSurface.addBilinear(
   innerProduct(Bs*dt(us), vs)
  +innerProduct((vel.expression()*grad(us)), vs)*Bs
  +innerProduct(divS(vel)*us, vs)*Bs
  +innerProduct(gradS(us), gradS(vs))*Ks*Bs
  // +innerProduct(Bb*ub-Bs*us, Bb*vb-Bs*vs)
  , interface
  , In
);



// R h3 = pow(meshSize,3);
// const & h_E(Parameter::meas);
bulkSurface.addFaceStabilization(
  innerProduct(1e-2*hi*jump(grad(ub)*n), jump(grad(vb)*n))
  ,Kh1
  ,In
);
bulkSurface.addBilinear(
  innerProduct(1e-2*jump(grad(us)*n), jump(grad(vs)*n))
  , Kh0
  , innerFacet
  , In
);

// impose initial condition
bulkSurface.addBilinear(
  innerProduct(ub,vb)*Bb
  , Kh1
);
bulkSurface.addBilinear(
  innerProduct(us, vs)*Bs
  , *interface(0)
  , In
  , 0
);

//       }
//       // std::cout << " Time Full assembly \t" << CPUtime() - tt0 << std::endl;
//       bulkSurface.addMatMul(data_uh);
//       // std::cout << " Time  A*u0 \t\t" << CPUtime() - tt0 << std::endl;
//
/*
CONSTRUCTION LINEAR PART
*/
// tt0 = CPUtime();
// impose initial condition
bulkSurface.addLinear(
  innerProduct(b0h.expression(), vb)*Bb
  , Kh1
);
bulkSurface.addLinear (
  innerProduct(s0h.expression(), vs)*Bs
  , *interface(0)
  , In , 0
);

// Add RHS on bulk
bulkSurface.addLinear(
        + innerProduct(fBh.expression(), Bb*vb)
        , Kh1
        , In
);

// Add RHS on surface
bulkSurface.addLinear(
        + innerProduct(fSh.expression(), Bs*vs)
        , interface
        , In
);

//       // std::cout << " Time assembly rhs \t" << CPUtime() - tt0 << std::endl;
//
//       /*
//       CONSTRUCTION NON LINEAR PART
//       */
//       bulkSurface.saveMatrix();
//       tt0 = CPUtime();
//       Expression2 wb(bh,0,op_id);
//       Expression2 ws(sh,0,op_id);
//
//       bulkSurface.addBilinear(
//         - innerProduct(wb*us , Bb*vb-Bs*vs)*Bbs
//         - innerProduct(ub*ws , Bb*vb-Bs*vs)*Bbs
//         , interface
//         , In
//       );
//
//       bulkSurface.addLinear(
//         - innerProduct(wb*ws , Bb*vb-Bs*vs)*Bbs
//         , interface
//         , In
//       );
//       // std::cout << " Time Full assembly Non linear part " << CPUtime() - tt0 << std::endl;
//     //
//     //   stokes.addLagrangeMultiplier(
//     //     innerProduct(1.,dp1), 0.
//     //     , In
//     //   );
//     //   stokes.addLagrangeMultiplier(
//     //     innerProduct(1.,dp1), 0.
//     //     , 0
//     //   );
//     //   // stokes.addLagrangeMultiplier(
//     //   //   innerProduct(1.,dp1), 0.
//     //   //   , 1
//     //   // );
//     //   // stokes.addLagrangeMultiplier(
//     //   //   innerProduct(1.,dp1), 0.
//     //   //   , 2
//     //   // );
//
      bulkSurface.solve();
//
//
//       Rn_ dwb(bulkSurface.rhs(SubArray(Wh.NbDoF()*In.NbDoF()   , 0)));
//       Rn_ dws(bulkSurface.rhs(SubArray(cutWh.NbDoF()*In.NbDoF(),idx_s0)));
//       double dist  = dwb.l1()/dwb.size();
//       double dists = dws.l1()/dws.size();
//       myCout << " [dwb , dws] = [ " << dist << " , " << dists << " ] " << std::endl;
//
Rn_ dw(bulkSurface.rhs_(SubArray(bulkSurface.get_nb_dof(), 0)));
data_uh = dw;
data_u0 = data_uh;
//
//       // Size has been changed by the lagrange multiplier
//       bulkSurface.rhs.resize(bulkSurface.nDoF); bulkSurface.rhs = 0.0;
//
//       iterNewton+=1;
//
//       if(iterNewton == 5 || max(dist, dists) < 1e-10 ) {
//
        bulkSurface.saveSolution(data_uh);
        bulkSurface.cleanMatrix();
//         break;
//       }
    // }
    // // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
    // {
    //   Fun2_h funX(Wh, fun_x);
    //   Fun2_h fun1(Wh, fun_1);
    //   R tt0 = CPUtime();
    //   double areaBubble   = integral(fun1, 0, 1) ;
    //   double centerOfMass = integral(funX, 1, 1) / areaBubble ;
    //
    //   double Pb = integralSurf(fun1, 1);
    //   double ra = sqrt(areaBubble / M_PI);
    //   double Pa = 2*M_PI*ra;
    //   double circularity = Pa / Pb;
    //   double riseVelocity = integral(uh, 1, 1) / areaBubble;
    //
    //   double q = integralSurf(us, 0);
    //
    //   std::cout << "\n Features of the drop " << std::endl;
    //   std::cout << " Time                   ->    " << GTime::current_time()+dT  << std::endl;
    //   std::cout << " Center Of Mass         ->    " << centerOfMass << std::endl;
    //   std::cout << " Circularity            ->    " << circularity << std::endl;
    //   std::cout << " Rise velocity          ->    " << riseVelocity << std::endl;
    //   std::cout << " Surfactant quantity    ->    " << q << std::endl;
    //
    //
    //   outputData << GTime::current_time()+dT << "\t"
    //              << centerOfMass << "\t"
    //              << circularity << "\t"
    //              << riseVelocity << "\t"
    //              << areaBubble <<  "\t"
    //              << q << std::endl;
    //
    // }

    // COMPUTE ERROR
    double errBulk =  L2normCut(b0h, fun_bulk2, GTime::current_time(), 0, 1);
    std::cout << " || uBulk-uBulkex||_2 = " << errBulk << std::endl;
    double errSurf =  L2normSurf(s0h, fun_surface, interface[0], GTime::current_time(), 0, 1);
    std::cout << " || uSurf-uSurfex||_2 = " << errSurf << std::endl;

    // errorsBulk.at(j) = errBulk;
    // errorsSurf.at(j) = errSurf;
    //
    //
    // // std::cout << " Set Velocity " << std::endl;
    // for(int i=0;i<nbTime;++i) {
    //   Fun2_h sol(Wh,In,data_uh);
    //   set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
    // }
    //
    // matlab::Export(bulkSurface.mat_, "matA.dat");
    // return 0;
    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
      std::cout << " Plotting " << std::endl;
      Fun_h errorBulk(Wh, fun_bulk, GTime::current_time());
      errorBulk.v -= b0h.v;
      errorBulk.v.map(fabs);
      Fun_h errorSurface(Sh, fun_surface, GTime::current_time());
      errorSurface.v -= s0h.v;
      errorSurface.v.map(fabs);

      Paraview<Mesh> writer(Kh1, "bulk_"+to_string(iterfig)+".vtk");
      writer.add(bh, "bulk", 0, 1);
      writer.add(errorBulk, "error", 0, 1);


      Paraview<Mesh> writerS(Kh0, "surface_"+to_string(iterfig)+".vtk");
      writerS.add(sh, "surface", 0, 1);
      writerS.add(errorSurface, "error", 0, 1);
      writerS.add(ls[0], "levelSet", 0, 1);

      // Paraview<Mesh> writerLS(Kh, "levelSet_"+to_string(iterfig)+".vtk");
      // writerLS.add(ls[0], "levelSet", 0, 1);


      iterfig++;
    }
    std::cout << " Timer iteration computation \t" << CPUtime() - tbegin_iter << std::endl;
    iter += 1;
  }
  // myCout << " -----------------------------------  " << std::endl;
  // myCout << " -----------------------------------  " << std::endl;
  // myCout << " END OF COMPUTATION " << std::endl;


  // Closing the opened files.
  // outputData.close();
}
