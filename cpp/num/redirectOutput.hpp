#ifndef REDIRECTOUTPUT_HPP_
#define REDIRECTOUTPUT_HPP_

#include <cstdlib>
#include <fstream>
#include <iostream>
// #include "../parallel/cfmpi.hpp"

// struct CoutFileAndScreen {
//   ofstream outFile;
//    ~CoutFileAndScreen(void){outFile.close();}

//   CoutFileAndScreen(std::string path2File) : outFile(path2File.c_str()) {
//     if(!MPIcf::IamMaster()) outFile.close();
//   }

//   CoutFileAndScreen& operator<< (ostream& (*pfun)(ostream&))
//    {
//      pfun(outFile);
//      pfun(cout);
//      return *this;
//    }
// };

// template <class T>
// CoutFileAndScreen& operator<< (CoutFileAndScreen& st, T val)
// {
//   if(st.outFile.is_open()) st.outFile << val;
//   cout << val;
//   return st;
// };

#endif
