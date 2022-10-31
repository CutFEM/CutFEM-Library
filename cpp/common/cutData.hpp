#ifndef _CUTDATA_HPP
#define _CUTDATA_HPP


struct CutData {
  typedef typename Mesh2::Element Element;

  bool cut = false;
  R2 vertices[2] = {R2(), R2()};
  byte sign[3];
  int edges[3] = {-1, -1, -1};

  R2 pointFromEdge(int i) const { assert(edges[i] != -1); return vertices[edges[i]];}

  CutData(const Element& T, double* s) {
    for(int i=0;i<3;++i) sign[i] = fsign(s[i]);
    if(s[0] != s[1] || s[1] != s[2]) {
      cut = true;
      int idv = 0;
      for(int i=0; i<3;++i) {
        const Ubyte v0 = Element::nvedge[i][0], v1= Element::nvedge[i][1];
        if(sign[v0] != sign[v1]) {
          assert(idv < 2);
          edges[i] = idv;
          const R t = -s[v0]/(s[v1]-s[v0]);
          vertices[idv] = (1-t) * ((R2) T[v0]) + t * ((R2) T[v1]) ;
          idv++;
        }
      }
    }
  }
  CutData(byte s[3]) {
    for(int i=0;i<3;++i) sign[i] = s[i];
  }
  CutData(byte s[3], R2 point[2], int arr_edge[3]) : CutData(s){
    if(s[0] != s[1] || s[1] != s[2]) {
      cut = true;
      vertices[0] = point[0];
      vertices[1] = point[1];
      edges[0] = arr_edge[0];
      edges[1] = arr_edge[1];
      edges[2] = arr_edge[2];
    }
  }
};






#endif
