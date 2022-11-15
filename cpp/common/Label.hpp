#ifndef LABEL_HPP
#define LABEL_HPP

#include <iostream>

class Label { // reference number for the physics
   friend inline std::ostream &operator<<(std::ostream &f, const Label &r) {
      f << r.lab;
      return f;
   }
   friend inline std::istream &operator>>(std::istream &f, Label &r) {
      f >> r.lab;
      return f;
   }

 public:
   int lab;
   Label(int r = 0) : lab(r) {}
   bool onGamma() const { return lab; }
};
#endif
