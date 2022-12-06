#ifndef REDIRECTOUTPUT_HPP_
#define REDIRECTOUTPUT_HPP_

#include <cstdlib>
#include <fstream>
#include <iostream>
// #include "../parallel/cfmpi.hpp"

struct CoutFileAndScreen {
   std::ofstream outFile;
   ~CoutFileAndScreen(void) { outFile.close(); }

   CoutFileAndScreen(std::string path2File) : outFile(path2File.c_str()) {
#ifdef USE_MPI
      if (!MPIcf::IamMaster())
         outFile.close();
#else
      outFile.close();
#endif
   }

   CoutFileAndScreen &operator<<(std::ostream &(*pfun)(std::ostream &)) {
      pfun(outFile);
      pfun(std::cout);
      return *this;
   }
};

template <class T> CoutFileAndScreen &operator<<(CoutFileAndScreen &st, T val) {
   if (st.outFile.is_open())
      st.outFile << val;
   std::cout << val;
   return st;
};

#endif
