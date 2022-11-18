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
