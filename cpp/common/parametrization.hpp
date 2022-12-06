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
#ifndef COMMON_PARAMETRIZATION_HPP
#define COMMON_PARAMETRIZATION_HPP

class CurveParametrization {

 public:
   virtual R2 evalutate(double t) const        = 0;
   virtual R2 evalutate(int i, double t) const = 0;
};

class LinearParametrization : public CurveParametrization {

 public:
   int nb_interval;

   struct LineSet {
      double a; // constant
      double b; // 1st order coefficient
      double x; // starting x value
   };

   LineSet *myLineSet;
   LinearParametrization(const std::vector<double> &x,
                         const std::vector<double> &y);
   void init(const std::vector<double> &x, const std::vector<double> &y);

   ~LinearParametrization() { delete[] myLineSet; }

 private:
   LinearParametrization(const LinearParametrization &);
   void operator=(const LinearParametrization &);
};

#endif
