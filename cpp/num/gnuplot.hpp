#ifndef GNUPLOT_HPP_
#define GNUPLOT_HPP_

#include <cstring>
#include <fstream>
#include "mapping.hpp"
#include "macroElement.hpp"

namespace gnuplot {

// void save(const Mesh2 & Th, std::string filename = "Th.dat") {

//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = Th[0].nv;
//   for(int k=0; k<Th.nt;++k) {
//     for(int i=0;i<nve;++i) {
//       plot << Th[k][i] << std::endl;
//     }
//     plot << Th[k][0] << std::endl;
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
// }

// void save(const MeshQuad2 & Th, std::string filename = "ThQ.dat") {

//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = Th[0].nv;
//   for(int k=0; k<Th.nt;++k) {
//     for(int i=0;i<nve;++i) {
//       plot << Th[k][i] << std::endl;
//     }
//     plot << Th[k][0] << std::endl;
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
// }
// void save(const Mesh3 & Th, std::string filename = "Th.dat") {

//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = Th[0].nv;
//   for(int k=0; k<Th.nt;++k) {
//     for(int i=0;i<nve;++i) {
//     plot << Th[k][i] << std::endl;
//     }
//     plot << Th[k][0] << std::endl;
//     plot << Th[k][2] << std::endl;
//     plot << Th[k][1] << std::endl;
//     plot << Th[k][3] << std::endl;
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
// }

template <typename Mesh>
void save(const TimeMacroElement<Mesh> &macro, std::string path = "./") {

   std::ofstream plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9,
       plot10;
   plot1.open("small1.dat", std::ofstream::out);
   plot2.open("small2.dat", std::ofstream::out);
   plot3.open("macroElement1.dat", std::ofstream::out);
   plot4.open("macroElement2.dat", std::ofstream::out);
   plot7.open("goodEdge1.dat", std::ofstream::out);
   plot8.open("goodEdge2.dat", std::ofstream::out);
   plot9.open("innerEdgeME1.dat", std::ofstream::out);
   plot10.open("innerEdgeME2.dat", std::ofstream::out);

   // const int nve = macro.Th[0].nbVertices();   // get number of vertices of
   // element
   const int nve = 3; // number of vertices of a triangle is three

   for (auto it = macro.small_element.begin(); it != macro.small_element.end();
        ++it) {

      int idxC = it->second.index; // index in cutSpace

      int domain = macro.Th.get_domain_element(idxC);
      // int idx = macro.Th.idxElementInBackMesh(idxC);

      if (domain == 0) {
         for (int i = 0; i < nve; ++i) {
            plot1 << macro.Th[i] << 1 << std::endl;
         }
         plot1 << macro.Th[0] << 1 << std::endl;
         plot1 << std::endl;
         plot1 << std::endl;
      } else {
         for (int i = 0; i < nve; ++i) {
            plot2 << macro.Th[i] << 2 << std::endl;
         }
         plot2 << macro.Th[0] << 2 << std::endl;
         plot2 << std::endl;
         plot2 << std::endl;
      }
   }
   plot1.close();
   plot2.close();
   int icolor = 0;
   for (auto it = macro.macro_element.begin(); it != macro.macro_element.end();
        ++it) {

      for (int i = 0; i < it->second.idx_element.size(); ++i) {

         int idxC   = it->second.idx_element[i];
         int domain = macro.Th.get_domain_element(idxC);
         // int idx = macro.Vh.idxElementInBackMesh(idxC);

         if (domain == 0) {

            for (int i = 0; i < nve; ++i) {
               plot3 << macro.Th[i] << icolor % 10 << std::endl;
            }
            plot3 << macro.Th[0] << icolor % 10 << std::endl;
            plot3 << std::endl;
            plot3 << std::endl;
         } else {
            for (int i = 0; i < nve; ++i) {
               plot4 << macro.Th[i] << icolor % 10 << std::endl;
            }
            plot4 << macro.Th[0] << icolor % 10 << std::endl;
            plot4 << std::endl;
            plot4 << std::endl;
         }
      }
      icolor += 3;
   }
   plot3.close();
   plot4.close();

   plot7.close();
   plot8.close();

   getchar();
   for (auto it = macro.macro_element.begin(); it != macro.macro_element.end();
        ++it) {

      for (int i = 0; i < it->second.inner_edge.size(); ++i) {

         int idxC   = it->second.inner_edge[i].first;
         int ie     = it->second.inner_edge[i].second;
         int domain = macro.Th.get_domain_element(idxC);
         int idx    = macro.Th.idxElementInBackMesh(idxC);

         int i0 = Mesh2::Element::nvedge[ie][0];
         int i1 = Mesh2::Element::nvedge[ie][1];
         R2 P   = 0.5 * (((R2)macro.Th[i0].H(idx)) + ((R2)macro.Th[i1].H(idx)));

         if (domain == 0) {
            plot9 << P << std::endl;

         } else {
            plot10 << P << std::endl;
         }
      }
   }
   plot9.close();
   plot10.close();
}

template <class Mesh> void save(const Mesh &mesh) {
   std::ofstream plot;
   plot.open("Th_nodes.dat", std::ofstream::out);
   for (int i = 0; i < mesh.nbVertices(); ++i) {
      plot << mesh(i).x << "\t" << mesh(i).y << std::endl;
   }
   plot.close();

   plot.open("Th_elements.dat", std::ofstream::out);
   for (int k = 0; k < mesh.nbElmts(); ++k) {
      const typename Mesh::Element &T(mesh[k]);
      for (int i = 0; i < 3; ++i)
         plot << mesh(T[i]) << "\t";
      plot << "0" << std::endl;
   }
   plot.close();

   plot.open("Th_edges.dat", std::ofstream::out);
   for (int k = 0; k < mesh.nbBrdElmts(); ++k) {
      const typename Mesh::BorderElement &T(mesh.be(k));
      for (int i = 0; i < 2; ++i)
         plot << mesh(T[i]) << "\t";
      for (int i = 0; i < 5; ++i)
         plot << "0 \t";
      plot << std::endl;
   }
   plot.close();
}

void save(const Interface<Mesh2> &Gh, std::string filename = "Gh.dat") {

   std::ofstream plot;
   plot.open(filename.c_str(), std::ofstream::out);
   for (int k = 0; k < Gh.nbElement(); ++k) {
      plot << Gh(k, 0) << std::endl;
      plot << Gh(k, 1) << std::endl;
      plot << std::endl;
      plot << std::endl;
   }

   plot.close();
}
void save(const Interface<Mesh2> &Gh, const Mapping2 &mapping,
          std::string filename = "Gh.dat") {

   std::ofstream plot;
   plot.open(filename.c_str(), std::ofstream::out);
   const int nve  = 2;
   const int accu = 10;

   for (int k = 0; k < Gh.nbElement(); ++k) {
      const int kb = Gh.idxElementOfFace(k); // on the back Mesh
      for (int i = 0; i < accu + 1; ++i) {
         R2 P(Gh(k, 0) * (double)i / accu + (1 - (double)i / accu) * Gh(k, 1));

         plot << mapping.map(kb, P) << std::endl;
      }
      plot << std::endl;
      plot << std::endl;
   }
   plot.close();
}

// void save(const Mesh2 & Th, const Fracture& fracture, std::string filename =
// "Th_fractured.dat") {
//
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = Th[0].nv;
//
//   Local_Partition local_partition;
//   for(int k=0; k<Th.nt;++k) {
//
//     // if(fracture.is_cut_element(k)) {
//       fracture.build_local_partition(k, local_partition);
//       // std::cout << " element " << k << " cut in " <<
//       local_partition.nb_element() << std::endl; for(int i=0;
//       i<local_partition.nb_element();++i){
//         Element2 K = local_partition.get_element(i);
//         for(int j=0;j<nve;++j) {
//           plot << K[j] << std::endl;
//         }
//         plot << K[0] << std::endl;
//         plot << std::endl;
//         plot << std::endl;
//       }
//       // for(int i=0;i<nve;++i) {
//       //   plot << Th[k][i] << std::endl;
//       // }
//       // plot << Th[k][0] << std::endl;
//       // plot << std::endl;
//       // plot << std::endl;
//       // getchar();
//     // }
//
//   }
//   plot.close();
// }
// void save(const Fracture& Gh, std::string filename = "fracture.dat") {
//
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = 2;
//   for(int k=0; k<Gh.nb_element();++k) {
//     for(int i=0;i<nve;++i) {
//       plot << Gh(k,i) << std::endl;
//     }
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
//
//   // plot.open("normal.dat", std::ofstream::out);
//   // for(int k=0; k<Gh.nbElement();++k) {
//   //   plot << 0.5*(Gh(k,0)+Gh(k,1)) << "\t" << 0.1*Gh.normal(k) <<
//   std::endl;
//   // }
//   plot.close();
//
// }
// void save(const Interface<Mesh2> & Gh, std::string filename = "Gh.dat") {
//
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = 2;
//   for(int k=0; k<Gh.nbElement();++k) {
//     for(int i=0;i<nve;++i) {
//   	plot << Gh(k,i) << std::endl;
//     }
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
//
//   plot.open("normal.dat", std::ofstream::out);
//   for(int k=0; k<Gh.nbElement();++k) {
//     plot << 0.5*(Gh(k,0)+Gh(k,1)) << "\t" << 0.1*Gh.normal(k) << std::endl;
//
//
//   }
//   plot.close();
//
// }

//
// void save(const Interface3 & Gh, std::string filename = "Gh.dat") {
//
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   const int nve = 3;
//   for(int k=0; k<Gh.nbElement();++k) {
//     for(int i=0;i<nve;++i) {
//   	plot << Gh(k,i) << std::endl;
//     }
//     plot << Gh(k,0) << std::endl;
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
// }

// void save(const Marker& marker, std::string filename = "marker.dat") {
//
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   for(int k=0; k<marker.size();++k) {
//     plot << marker.get_marker(k) << std::endl;
//   }
//   plot.close();
//
//   // plot.open("node.dat", std::ofstream::out);
//   // const int nve = 2;
//   // for(int k=0; k<marker.nbElement();++k) {
//   //   for(int i=0;i<nve;++i) {
//   // 	plot << marker(k,i) << std::endl;
//   //   }
//   //   plot << std::endl;
//   //   plot << std::endl;
//   // }
//   // plot.close();
//
// }

// void saveNormal(const Marker & marker, std::string filename = "normal.dat") {
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   // for(int k=0; k<marker.edges_node.size();++k) {
//   //   plot << marker.edges_node[k] << std::endl;
//   // }
//   // plot << marker.edges_node[0] << std::endl;
//   for(int k=0; k<marker.faces_.size();++k) {
//     const typename Marker::Face& face = marker[k];  // the face
//     R2 normal(marker.normal(k));
//     plot << (marker(k,0)+marker(k,1))*0.5 << "\t" << 0.1*normal << std::endl;
//   }
//   plot.close();
//
// }

//
// void save(const MacroElement & macro, const Extension& extension) {
//
//   std::ofstream plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,
//   plot9, plot10, plot11, plot12; plot1.open("small1.dat",
//   std::ofstream::out); plot2.open("small2.dat", std::ofstream::out);
//   plot3.open("macroElement1.dat", std::ofstream::out);
//   plot4.open("macroElement2.dat", std::ofstream::out);
//   plot5.open("extensionEdge1.dat", std::ofstream::out);
//   plot6.open("extensionEdge2.dat", std::ofstream::out);
//   plot7.open("goodEdge1.dat", std::ofstream::out);
//   plot8.open("goodEdge2.dat", std::ofstream::out);
//   plot11.open("exhaustEdge1.dat", std::ofstream::out);
//   plot12.open("exhaustEdge2.dat", std::ofstream::out);
//   plot9.open("innerEdgeME1.dat", std::ofstream::out);
//   plot10.open("innerEdgeME2.dat", std::ofstream::out);
//
//   const int nve = macro.Vh.Th[0].nv;
//   for(auto it=macro.small_element.begin();
//   it!=macro.small_element.end();++it) {
//
//     int idxC = it->second.index;  // index in cutSpace
//     int domain = macro.Vh.whichDomain(idxC);
//     int idx = macro.Vh.idxElementInBackMesh(idxC);
//
//     if(domain == 0){
//       for(int i=0;i<nve;++i) {
//         plot1 << macro.Vh.Th[idx][i] << 1 << std::endl;
//       }
//       plot1 << macro.Vh.Th[idx][0] << 1 << std::endl;
//       plot1 << std::endl;
//       plot1 << std::endl;
//     }
//     else {
//       for(int i=0;i<nve;++i) {
//         plot2 << macro.Vh.Th[idx][i] << 2 << std::endl;
//       }
//       plot2 << macro.Vh.Th[idx][0] << 2 << std::endl;
//       plot2 << std::endl;
//       plot2 << std::endl;
//     }
//
//   }
//   plot1.close();
//   plot2.close();
//   int icolor = 0;
//   for(auto it=macro.macro_element.begin();
//   it!=macro.macro_element.end();++it) {
//
//     // std::cout << " Macro \t" << it->second.idx_root_element << std::endl;
//     for(int i=0;i<it->second.idx_element.size();++i) {
//
//       int idxC = it->second.idx_element[i];
//       int domain = macro.Vh.whichDomain(idxC);
//       int idx = macro.Vh.idxElementInBackMesh(idxC);
//
//       if(domain == 0){
//
//         for(int i=0;i<nve;++i) {
//           plot3 << (R2)macro.Vh.Th[idx][i] << "\t" << icolor%10 << std::endl;
//         }
//         plot3 << (R2)macro.Vh.Th[idx][0] << "\t" << icolor%10 << std::endl;
//         plot3 << std::endl;
//         plot3 << std::endl;
//       }
//       else {
//         // std::cout << idxC << std::endl;
//         for(int i=0;i<nve;++i) {
//           plot4 << (R2)macro.Vh.Th[idx][i] << "\t" << icolor%10 << std::endl;
//         }
//         plot4 << (R2)macro.Vh.Th[idx][0] << "\t" << icolor%10 << std::endl;
//         plot4 << std::endl;
//         plot4 << std::endl;
//       }
//     }
//     icolor+= 3;
//   }
//   plot3.close();
//   plot4.close();
//   for(auto it=extension.element_edge_handle.begin();
//   it!=extension.element_edge_handle.end();++it) {
//
//     int idxK = it->first.first;
//     int ie = it->first.second;
//     int domain = macro.Vh.whichDomain(idxK);
//     int idx = macro.Vh.idxElementInBackMesh(idxK);
//     int handle = it->second;
//     int i0 = Mesh2::Element::nvedge[ie][0];
//     int i1 = Mesh2::Element::nvedge[ie][1];
//     R2 P = 0.5*(((R2) macro.Vh.Th[idx][i0]) + ((R2) macro.Vh.Th[idx][i1]));
//     if(handle == 1) {
//       if(domain == 0){
//         plot5 << P << std::endl;
//       }
//       else {
//         plot6 << P << std::endl;
//       }
//     }
//     else if(handle == 0){
//       if(domain == 0){
//         plot7 << P << std::endl;
//       }
//       else {
//         plot8 << P << std::endl;
//       }
//     }
//     else{
//       if(domain == 0){
//         plot11 << P << std::endl;
//       }
//       else {
//         plot12 << P << std::endl;
//       }
//
//     }
//   }
//   plot5.close();
//   plot6.close();
//   plot7.close();
//   plot8.close();
//   plot11.close();
//   plot12.close();
//   for(auto it=macro.macro_element.begin();
//   it!=macro.macro_element.end();++it) {
//
//     for(int i=0;i<it->second.inner_edge.size();++i) {
//
//       int idxC = it->second.inner_edge[i].first;
//       int ie = it->second.inner_edge[i].second;
//       int domain = macro.Vh.whichDomain(idxC);
//       int idx = macro.Vh.idxElementInBackMesh(idxC);
//
//       int i0 = Mesh2::Element::nvedge[ie][0];
//       int i1 = Mesh2::Element::nvedge[ie][1];
//       R2 P = 0.5*(((R2) macro.Vh.Th[idx][i0]) + ((R2) macro.Vh.Th[idx][i1]));
//
//       if(domain == 0){
//         plot9 << P << std::endl;
//
//       }
//       else {
//         plot10 << P << std::endl;
//       }
//     }
//   }
//   plot9.close();
//   plot10.close();
//
// }
//
// void save(const MacroElementSurface & macro) {
//   std::ofstream plot;
//   plot.open("small.dat", std::ofstream::out);
//
//   const Mesh2& Th(*macro.interface.backMesh);
//   const int nve = Th[0].nv;
//   for(auto it=macro.small_element.begin();
//   it!=macro.small_element.end();++it) {
//
//     int idxC = it->first;  // index in cutMesh
//     int idx = macro.interface.idxElementOfFace(idxC);
//
//     for(int i=0;i<nve;++i) {
//       plot << Th[idx][i] << 1 << std::endl;
//     }
//     plot << Th[idx][0] << 1 << std::endl;
//     plot << std::endl;
//     plot << std::endl;
//   }
//   plot.close();
//
//
//   plot.open("macroElement.dat", std::ofstream::out);
//   int icolor = 0;
//   for(auto it=macro.macro_element.begin();
//   it!=macro.macro_element.end();++it) {
//     for(int i=0;i<it->second.idx_element.size();++i) {
//
//       int idxC = it->second.idx_element[i];
//       int idx = macro.interface.idxElementOfFace(idxC);
//
//       for(int i=0;i<nve;++i) {
//         plot << Th[idx][i] << icolor%10 << std::endl;
//       }
//       plot << Th[idx][0] << icolor%10 << std::endl;
//       plot << std::endl;
//       plot << std::endl;
//     }
//     icolor+= 3;
//   }
//   plot.close();
//
//   plot.open("innerEdge.dat", std::ofstream::out);
//   for(auto it=macro.macro_element.begin();
//   it!=macro.macro_element.end();++it) {
//
//     for(int i=0;i<it->second.inner_edge.size();++i) {
//
//       int iface = it->second.inner_edge[i].first;
//       int ie = it->second.inner_edge[i].second;
//       int idx = macro.interface.idxElementOfFace(iface);
//
//       int i0 = Mesh2::Element::nvedge[ie][0];
//       int i1 = Mesh2::Element::nvedge[ie][1];
//       R2 P = 0.5*(((R2) Th[idx][i0]) + ((R2) Th[idx][i1]));
//
//       plot << P << std::endl;
//
//     }
//   }
//   plot.close();
//
//
// }
//

} // namespace gnuplot

#endif
