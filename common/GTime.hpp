#ifndef GTIME_HPP
#define GTIME_HPP

#include <iostream>
#include <iomanip>
#include <cassert>

class GTime {
public :
    
  // static double current_time;
  static double time_step;
  static int current_iteration;
  static int total_number_iteration;
  static int iter_step_quadrature;
  static double t0;
  
  static void info() {
    std::cout << " INFO OF THE TIME STEPPING " << std::endl;
    std::cout << std::setw(30) << std::left <<" current time "
	      << std::setw(10) << " " << current_time() << std::endl;
    std::cout << std::setw(30) << std::left <<" current iteration "
	      << std::setw(10) << " " << current_iteration << std::endl;
    std::cout << std::setw(30) << std::left <<" time step "
	      << std::setw(10) << " " << time_step << std::endl;
    std::cout << std::setw(30) << std::left <<" total number of iteration "
	      << std::setw(10) << " " << total_number_iteration << std::endl;
  }

  static double final_time() {return t0 + time_step*total_number_iteration;}
  static double current_time() {return t0 + time_step*current_iteration;}
  
};



// double GTime::current_time = 0;
// double GTime::time_step = 0;
// int GTime::total_number_iteration = 0;
// int GTime::current_iteration = 0;
// int GTime::iter_step_quadrature = 0;
// double GTime::t0 = 0;
#endif
