#ifndef GNUPLOT_HPP_
#define GNUPLOT_HPP_

#include <cstring>
#include <fstream>

namespace gnuplot {



  void save(const Mesh2 & Th, std::string filename = "Th.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    const int nve = Th[0].nv;
    for(int k=0; k<Th.nt;++k) {
      for(int i=0;i<nve;++i) {
        plot << Th[k][i] << std::endl;
      }
      plot << Th[k][0] << std::endl;
      plot << std::endl;
      plot << std::endl;
    }
    plot.close();
  }

  void save(const CutFESpace2 & Vh, int domain, std::string filename = "CutTh.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    const int nve = Vh.Th[0].nv;
    for(int k=0; k<Vh.NbElement(domain);++k) {
      int kk = Vh.subDomain(domain)->getTriLocToGlob(k);
      for(int i=0;i<nve;++i) {
        plot << Vh.Th[kk][i] << std::endl;
      }
      plot << Vh.Th[kk][0] << std::endl;
      plot << std::endl;
      plot << std::endl;
    }
    plot.close();
  }


  void save(const Mesh3 & Th, std::string filename = "Th.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    const int nve = Th[0].nv;
    for(int k=0; k<Th.nt;++k) {
      for(int i=0;i<nve;++i) {
    	plot << Th[k][i] << std::endl;
      }
      plot << Th[k][0] << std::endl;
      plot << Th[k][2] << std::endl;
      plot << Th[k][1] << std::endl;
      plot << Th[k][3] << std::endl;
      plot << std::endl;
      plot << std::endl;
    }
    plot.close();
  }


  void save(const Interface2 & Gh, std::string filename = "Gh.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    const int nve = 2;
    for(int k=0; k<Gh.nbElement();++k) {
      for(int i=0;i<nve;++i) {
    	plot << Gh(k,i) << std::endl;
      }
      plot << std::endl;
      plot << std::endl;
    }
    plot.close();

    plot.open("normal.dat", std::ofstream::out);
    for(int k=0; k<Gh.nbElement();++k) {
      plot << 0.5*(Gh(k,0)+Gh(k,1)) << "\t" << 0.1*Gh.normal(k) << std::endl;


    }
    plot.close();

  }

  void save(const Interface3 & Gh, std::string filename = "Gh.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    const int nve = 3;
    for(int k=0; k<Gh.nbElement();++k) {
      for(int i=0;i<nve;++i) {
    	plot << Gh(k,i) << std::endl;
      }
      plot << Gh(k,0) << std::endl;
      plot << std::endl;
      plot << std::endl;
    }
    plot.close();
  }



  void saveMarker(const Marker & marker, std::string filename = "marker.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    for(int k=0; k<marker.markers.size();++k) {
      plot << marker.markers[k] << std::endl;
    }
    // plot << marker.markers[0] << std::endl;
    plot.close();

    plot.open("node.dat", std::ofstream::out);
    const int nve = 2;
    for(int k=0; k<marker.nbElement();++k) {
      for(int i=0;i<nve;++i) {
    	plot << marker(k,i) << std::endl;
      }
      plot << std::endl;
      plot << std::endl;
    }
    plot.close();

  }

  void saveNormal(const Marker & marker, std::string filename = "normal.dat") {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    // for(int k=0; k<marker.edges_node.size();++k) {
    //   plot << marker.edges_node[k] << std::endl;
    // }
    // plot << marker.edges_node[0] << std::endl;
    for(int k=0; k<marker.faces_.size();++k) {
      const typename Marker::Face& face = marker[k];  // the face
      R2 normal(marker.normal(k));
      plot << (marker(k,0)+marker(k,1))*0.5 << "\t" << 0.1*normal << std::endl;
    }
    plot.close();

  }

  // void exportInterface(const Interface2 & Gh, const Mapping2& mapping,
  // 		       std::string filename = "Gh.dat") {

  //   std::ofstream plot;
  //   plot.open(filename.c_str(), std::ofstream::out);
  //   const int nve = 2;
  //   const int accu = 10;

  //   for(int k=0; k<Gh.nbElement();++k) {

  //     const int kb = Gh.idxElementOfFace(k);                // on the back Mesh
  //     const int kl = mapping.Vh.idxInLocMesh(kb);           // on mapping Mesh
  //     const typename Mapping2::FElement & FK(mapping.Vh[k]);

  //     for(int i=0;i<accu+1;++i) {
  // 	R2 P(Gh(k,0)*(double)i/accu + (1 - (double)i/accu)*Gh(k,1));
  // 	  plot << mapping.deformNode(FK, P) << std::endl;

  //     }
  //     plot << std::endl;
  //     plot << std::endl;
  //   }
  //   plot.close();


  // }




  void save(const MacroElement & macro) {

    std::ofstream plot1, plot2;
    plot1.open("small1.dat", std::ofstream::out);
    plot2.open("small2.dat", std::ofstream::out);

    const int nve = macro.Vh.Th[0].nv;
    for(int k=0; k<macro.idx_small_K.size();++k) {

      int idxC = macro.idx_small_K[k].first;  // index in cutSpace
      int domain = macro.Vh.whichDomain(idxC);
      int idx = macro.Vh.idxElementInBackMesh(idxC);

      if(domain == 0){
        for(int i=0;i<nve;++i) {
          plot1 << macro.Vh.Th[idx][i] << 1 << std::endl;
        }
        plot1 << macro.Vh.Th[idx][0] << 1 << std::endl;
        plot1 << std::endl;
        plot1 << std::endl;
      }
      else {
        for(int i=0;i<nve;++i) {
          plot2 << macro.Vh.Th[idx][i] << 2 << std::endl;
        }
        plot2 << macro.Vh.Th[idx][0] << 2 << std::endl;
        plot2 << std::endl;
        plot2 << std::endl;
      }

    }
    plot1.close();
    plot2.close();
  }




}


#endif
