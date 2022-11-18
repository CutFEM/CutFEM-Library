#ifndef COMMON_GLOBAL_HPP
#define COMMON_GLOBAL_HPP

#include <numeric>
namespace globalVariable {
extern int verbose;
extern double UnSetMesure;
// const double DoubleEpsC = 1.0e-9; // numeric_limits<double>::epsilon();
const double Epsilon = 10 * std::numeric_limits<double>::epsilon();
// const double pi = 3.14159265359;
}; // namespace globalVariable

#endif
