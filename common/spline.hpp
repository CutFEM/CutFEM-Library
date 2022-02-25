#ifndef _SPLINE_HPP_
#define _SPLINE_HPP_
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "R2.hpp"
using namespace std;

class Spline2D;

//  A 2D cubic spline
class Spline {
public:
  Spline(){};
  Spline( const vector<double> & x, const vector<double> & y );
private:
  // The fitted points
  vector<double> X;
  vector<double> Y;

  /// The coefficients of the spline curve between two points
    struct SplineSet
    {
      double a;   // constant
      double b;   // 1st order coefficient
      double c;   // 2nd order coefficient
      double d;   // 3rd order coefficient
      double x;   // starting x value
      void print() const {
        std::cout << " x_0 = " << x << std::endl;
        std::cout << "p(x) = " << a << " + "
                               << b << "x + "
                               << c << "x^2 + "
                               << d << "x^3 " << std::endl;
      }
    };

  // The coefficients of the spline curves between all points
  SplineSet * mySplineSet;

  void init( const vector<double> & x, const vector<double> & y );

public:
  /** Get the Y value of the spline curves for a particular X
    input x
    return the y value
  */
  double evaluate(const double) const;
  double operator()(const double a) const {return evaluate(a);}
  /* Compute the y value of the spline derivative for a particular x
    input x
    return the value y of the diff in x
   */
  double diff(const double) const;
  double ddiff(const double) const;
  /* Write in a file the spline to display it with gnuplot
    input N: number of points between two points of the set
    input filename: the files in which the solution is written
    output err: error in the computation or not
   */
  int gnuplot(int , string );
  friend  void display_spline2d(const Spline &, const Spline &, string, int);
  friend class Spline2D;
   ~Spline(){ delete [] mySplineSet;}
private:
  Spline(const Spline &);
  void operator=(const Spline &);

};

class Spline2D {


  Spline splineX;
  Spline splineY;

  int nb_node;

public:
  Spline2D( const vector<double> & t, const vector<double> & x, const vector<double> & y )
  : splineX(t, x), splineY(t,y) {
    nb_node = t.size();
    assert((x.size() == nb_node) && (y.size() == nb_node));
  }
  Spline2D( const vector<double> & t, const vector<R2> & x );

public:
  R2 evaluate(double t) const ;
  int gnuplot(int , string );

  // ~Spline2D(){delete []splineX; delete splineY;}
private:
  Spline2D(const Spline2D &);
  void operator=(const Spline2D &);
};

void display_spline2d(const Spline &, const Spline &, string, int=10);



#endif
