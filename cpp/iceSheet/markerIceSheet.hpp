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
#ifndef _MARKER_ICESHEET_HPP
#define _MARKER_ICESHEET_HPP

#include "../common/marker.hpp"

class MarkerIceSheet : public Marker {

 public:
   R2 groundingNode;
   int idxGroundingNode;

   MarkerIceSheet(const Mesh2 &Thh);

   void setLabelBedRock();
   void setBoundaryFace();
   void findGroundingPoint();
};

#endif
