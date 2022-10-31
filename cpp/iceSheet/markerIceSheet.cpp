#include "markerIceSheet.hpp"
#include "../common/dataSara.hpp"



MarkerIceSheet::MarkerIceSheet(const Mesh2& Thh) : Marker(Thh) {


  DataSara data("../../iceSheet/");

  add("../../iceSheet/surf.dat", 1);   // label = 1
  add("../../iceSheet/bed.dat");
  setLabelBedRock();
  setDomainSign(data.ls_sign);
  setBoundaryFace();

}



void MarkerIceSheet::setLabelBedRock() {

  findGroundingPoint();

  for(int i=nstartFace[1]; i<nendFace[1]; ++i) {

    FaceMarker& face(faces[i]);
    int i1 = face[0];
    int i2 = face[1];

    if(i1 < idxGroundingNode) {
      face.lab = 2;
    }
    else {
      face.lab = 3;
    }
  }

}


// Node with minimum y-value
void MarkerIceSheet::findGroundingPoint() {

  groundingNode = edges_node[nstart[1]];
  for(int i=nstartMarker[1]; i<nendMarker[1];++i) {
    R2 P = edges_node[i];
    if(P.y <  groundingNode.y) {
      groundingNode = P;
      idxGroundingNode = i;
    }

  }
}


void MarkerIceSheet::setBoundaryFace() {

  for(int k=0;k<Th.nt;++k) {

    if(isNotCut(k) && domain_sign(Th(k,0)) > 0) continue; // only Omega2

    const Element& K(Th[k]);

    for(int i=0;i<Element::ne;++i) {
      int jn = i;
      int k_next = Th.ElementAdj(k, jn);

      if(k_next != -1) continue;  // only want border element

      int ie[2] = {Element::nvedge[i][1] ,Element::nvedge[i][0]};

      if(domain_sign(Th(k,ie[0])) > 0 && domain_sign(Th(k,ie[1])) > 0) continue; // not cut && not in domain

      if(isCut(k)) {
        // add the Cut Face
        int iface = element_seen(k);
        FaceMarker& face(faces[iface]);
        int idx0 = -1;
        // find the point on edge
        for(int j=0;j<2;++j) {
          R2 P = (edges_node[face[j]]);
          R2 ed1(P,K[ie[0]]);
          R2 ed2(P,K[ie[1]]);
          if( fabs((ed1,ed2) + ed1.norm()*ed2.norm()) < 100*Epsilon)
          idx0 = j;
        }

        int idx1 = -1;
        for(int j=0;j<2;++j) {
          if(domain_sign(Th(k,ie[j])) < 0)
          idx1 = ie[j];
        }

        assert(idx0!= -1);
        assert(idx1!= -1);

        R2 A0 = edges_node[face[idx0]];
        R2 A1 = K[idx1];

        if(((A1-A0).perp(),K.N(i)) < 0) {
          R2 As = A0;
          A0 = A1;
          A1 = As;
        }

        edges_node.push_back(A0);
        edges_node.push_back(A1);
        int i1 = edges_node.size()-2;
        int i2 = edges_node.size()-1;
        int labb = (A0.x < 1e-5)? 4: 5;
        faces.push_back(FaceMarker(*this, k, i1,i2, labb));
        cut_element.push_back(k);
      }
      else {
        // add the edge
        edges_node.push_back(K[ie[0]]);
        edges_node.push_back(K[ie[1]]);
        int i1 = edges_node.size()-2;
        int i2 = edges_node.size()-1;
        int labb = (K[ie[0]].x < 1e-5)? 4: 5;
        faces.push_back(FaceMarker(*this, k, i1,i2, labb));
        cut_element.push_back(k);
      }
    }
  }
}
