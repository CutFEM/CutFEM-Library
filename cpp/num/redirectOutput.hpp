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
        if (!MPIcf::IamMaster())
            outFile.close();
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
