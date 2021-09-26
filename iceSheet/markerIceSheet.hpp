#ifndef _MARKER_ICESHEET_HPP
#define _MARKER_ICESHEET_HPP

#include "../common/marker.hpp"


class MarkerIceSheet : public Marker {



public:
  R2 groundingNode;
  int idxGroundingNode;
  
  MarkerIceSheet(const Mesh2& Thh); 


  void setLabelBedRock();
  void setBoundaryFace();
  void findGroundingPoint();







};




#endif
