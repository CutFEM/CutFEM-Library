#include "matlab.hpp"

void matlab::Export(std::map<std::pair<int,int>,R> & dF, std::string filename) {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    plot << std::setprecision(18);
    for (std::map<std::pair<int,int>,R>::const_iterator q=dF.begin();
	 q != dF.end(); ++q)
      {
	plot << q->first.first  + 1 << "\t" << q->first.second + 1
	     << "\t" << q->second << std::endl;
      }
    plot.close();
  }
