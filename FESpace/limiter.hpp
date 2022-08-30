#ifndef CUTFEM_LIMITER_HPP
#define CUTFEM_LIMITER_HPP


#include "FESpace.hpp"
#include "macroElement.hpp"




namespace limiter {

  // check mean value to satisfy maximuim principle
  template<typename Mesh>
  void check_maximum_principle(std::map<int, double>& u_mean, double min_u, double max_u)  {

    for (auto& p : u_mean) {
      double val = p.second;
      if( min_u > val || max_u < val) {
        std::cout << "element \t" << p.first << "\t" << min_u
        << "\t" << val << "\t" << max_u << std::endl;
      }
    }
  }


  namespace CutFEM {

    // -------------------------------------------------------------------------
    // COMPUTE THE MIN AND THE MAX VALUE OF A FE FUNCTION ON A CUT MESH
    // DEPENDING OF THE POLYNOMIAL ORDER.
    // BY DEFAULT IT USES LINEAR ELEMENTS

    template<typename Mesh>
    void findMinAndMaxValue_P0(const FunFEM<Mesh>& uh, double& min_val, double& max_val) {
      typedef typename Mesh::Rd Rd;
      typedef typename Mesh::Element Element;
      typedef typename GFESpace<Mesh>::FElement FElement;
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const ActiveMesh<Mesh>& Kh(Wh.get_mesh());

      double min_val_loc=1e300;
      double max_val_loc = -1e300;
      int k_min, k_max;

      for (int k=Wh.first_element();k<Wh.last_element();
             k+= Wh.next_element()) {
        const FElement& FK(Wh[k]);

        Rd mip = FK.T.barycenter();
        double val = uh.eval(k, mip);
        if(val > max_val_loc){
          max_val_loc = val; k_max = k;
        }
        if(val < min_val_loc) {
          min_val_loc = val; k_min=k;
        }
      }
      MPIcf::AllReduce(min_val_loc, min_val, MPI_MIN);
      MPIcf::AllReduce(max_val_loc, max_val, MPI_MAX);
    }

    template<typename Mesh>
    void findMinAndMaxValue_P1(const FunFEM<Mesh>& uh, double& min_val, double& max_val) {
      typedef typename Mesh::Rd Rd;
      typedef typename Mesh::Element Element;
      typedef typename GFESpace<Mesh>::FElement FElement;
      int nv_loc = Mesh::Rd::d +1;
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const ActiveMesh<Mesh>& Kh(Wh.get_mesh());

      // min_val=1e300;
      // max_val = -1e300;
      double min_val_loc=1e300;
      double max_val_loc = -1e300;
      // Apply limiter in all element
      int k_min, k_max;
      // for(int k=Wh.first_element(); k<Wh.last_element();++k) {
      for (int k=Wh.first_element();k<Wh.last_element();
             k+= Wh.next_element()) {
        const FElement& FK(Wh[k]);
        const Cut_Part<Element> cutK(Kh.get_cut_part(k,0));


        for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
          // get the local min and max on K
          for(int ipq = 0; ipq < nv_loc; ++ipq)  {
            Rd mip = cutK.get_vertex(it, ipq);
            // double val = uh.eval(k, mip);
            double val = uh.eval(k, mip);
            if(val > max_val_loc){
              max_val_loc = val; k_max = k;
            }
            if(val < min_val_loc) {
              min_val_loc = val; k_min=k;
            }
            // max_val  = max(max_val, val);
            // min_val  = min(min_val, val);
          }
        }

      }
      MPIcf::AllReduce(min_val_loc, min_val, MPI_MIN);
      MPIcf::AllReduce(max_val_loc, max_val, MPI_MAX);
      // std::cout << " minimum on element " << k_min << std::endl;
      // std::cout << " maximum on element " << k_max << std::endl;
    }

    template<typename Mesh>
    void findMinAndMaxValue(const FunFEM<Mesh>& uh, double& min_val, double& max_val) {
      const int fctOrder = uh.getPolynomialOrder();
      if(fctOrder == 0) {
        findMinAndMaxValue_P0(uh, min_val, max_val);
      }
      else if (fctOrder == 1 ){
        findMinAndMaxValue_P1(uh, min_val, max_val);
      }
      else {
        findMinAndMaxValue_P1(uh, min_val, max_val);
        std::cout << " Min and max compute using linear elements. Not exact for this polynomial order" << std::endl;
      }

    }


    // -----------------------------------------------------------------------------


    // Only for discontinuous Lagrange
    template<typename Mesh>
    void extendToMacro_P0(const FunFEM<Mesh>& uh, Rn& u_new, std::map<int, double>& u_mean, const MacroElement<Mesh>& macro) {
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename Mesh::Rd Rd;
      u_new = uh.v;
      const GFESpace<Mesh>& Wh(*uh.Vh);
      for (auto& p : macro.macro_element) {

        // compute the value of the dof of the elements in the
        // macro element.
        // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K \cap \Omega|}
        const MElement& MK(p.second);
        int n_element = MK.size();

        // loop over the elements that has to be changed
        for(auto s : MK.idx_element ) {
          const FElement& FK(Wh[s]);
          assert(FK.NbDoF() == FK.NbNode());
          double area_total = MK.area_total_;

          // loop over the dof (in that case just node evaluation)
          for(int df=0; df<FK.NbNode(); ++df) {
            Rd P = FK.Pt(df);
            int df_glb = FK.loc2glb(df);
            double val = 0;

            for(auto k : MK.idx_element ) {
              double s = macro.get_area(k);
              val += s * uh.eval(k, P);
            }
            val = val / area_total;
            u_new(df_glb) = val;
          }
        }

      }
    }

    // Only for discontinuous Lagrange
    template<typename Mesh>
    void extendToMacro_P1(const FunFEM<Mesh>& uh, Rn& u_new, std::map<int, double>& u_mean, const MacroElement<Mesh>& macro) {
      typedef typename Mesh::Element Element;
      bool mean_good = true;
      double mm = 0;
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename FElement::QF QF;
      typedef typename Mesh::Rd Rd;

      const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
      const GFESpace<Mesh>& Wh(*uh.Vh);
      FunFEM<Mesh> fun_u_new(Wh, u_new);

      for (auto& p : macro.macro_element) {

        // compute the value of the dof of the elements in the
        // macro element.
        // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K \cap \Omega|}
        const MElement& MK(p.second);
        int n_element = MK.size();
        double area_total = MK.area_total_;
        int idx_root = MK.idx_root_element;

        // loop over the elements that has to be changed
        for(auto s : MK.idx_element ) {
          const FElement& FK(Wh[s]);

          // loop over the dof (in that case just node evaluation)
          for(int df=0; df<FK.NbDoF() ; ++df) {
            Rd P = FK.Pt(df);
            int df_glb = FK.loc2glb(df);
            double val = 0;
            for(auto k : MK.idx_element ) {
              double ss = macro.get_area(k);
              val += ss * uh.eval(k, P);
            }
            val = val / area_total;
            u_new(df_glb) = val;

          }
        }

        double C0 = 0;
        double mean_val = 0.;
        for(auto k : MK.idx_element ) {
          const FElement& FK(Wh[k]);
          const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));

          for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
            double meas = cutK.measure(it);

            for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
              typename QF::QuadraturePoint ip(qf[ipq]); // integration point
              Rd mip = cutK.mapToPhysicalElement(it, ip);
              double Cint = meas * ip.getWeight();
              double val_old = uh.eval(k, mip);
              mean_val += Cint*val_old;
              C0 += Cint*(val_old - fun_u_new.eval(k,mip));
              // if(val_old < 0 ){
              //
              //   std::cout << "-------- \t" << val_old << std::endl;
              //   for(int df=0; df<FK.NbDoF() ; ++df) {
              //     std::cout << uh.v[df] << std::endl;
              //   }
              //   getchar();
              // }
            }
          }
        }
        C0 = C0 / area_total;
        mean_val = mean_val / area_total;
        // if(mean_val < 0 || mean_val > 1.) {
        //     // std::cout << mean_val << std::endl;
        //     // getchar();
        //   mm = mean_val;
        //   mean_good = false;
        // }
        u_mean[idx_root] = mean_val;

        for(auto k : MK.idx_element ) {
          const FElement& FK(Wh[k]);
          for(int df=0; df<FK.NbDoF() ; ++df) {
            int df_glb = FK.loc2glb(df);
            u_new(df_glb) += C0;
          }
        }



        // double mean1=0., mean2=0.;
        // for(auto k : MK.idx_element ) {
        //   const FElement& FK(Wh[k]);
        //   const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
        //
        //   for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
        //     double meas = cutK.measure(it);
        //
        //     for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
        //       typename QF::QuadraturePoint ip(qf[ipq]); // integration point
        //       Rd mip = cutK.mapToPhysicalElement(it, ip);
        //       double Cint = meas * ip.getWeight();
        //       mean1 += Cint*uh.eval(k, mip) / area_total;
        //       mean2 += Cint*fun_u_new.eval(k, mip) / area_total;
        //     }
        //   }
        // }
        // std::cout << mean1 << "\t" << mean2 << "\t" << u_mean[idx_root] << std::endl;
        // assert(fabs(mean1-mean2)< 1e-15);
        // assert(fabs(u_mean[idx_root]-mean2)< 1e-15);
      }
      if(!mean_good) {
        std::cout << " Mean value not good \t" << mm  << std::endl;
      }
    }

    template<typename Mesh>
    void extendToMacro(const FunFEM<Mesh>& uh, Rn& u_new, std::map<int, double>& u_mean, const MacroElement<Mesh>& macro) {
      const BasisFctType basisFctType = uh.getBasisFctType();
      if(basisFctType == BasisFctType::P0) {
        extendToMacro_P0(uh, u_new, u_mean, macro);
      }
      else if (basisFctType == BasisFctType::P1dc ){
        extendToMacro_P1(uh, u_new, u_mean, macro);
      }
      else {
        std::cout << " Extension to macro element not implemented for those element" << std::endl;
      }
    }
    template<typename Mesh>
    void extendToMacro(const FunFEM<Mesh>& uh, Rn& u_new, const MacroElement<Mesh>& macro) {
      const BasisFctType basisFctType = uh.getBasisFctType();
      std::map<int, double> u_mean;
      if(basisFctType == BasisFctType::P0) {
        extendToMacro_P0(uh, u_new, u_mean, macro);
      }
      else if (basisFctType == BasisFctType::P1dc ){
        extendToMacro_P1(uh, u_new, u_mean, macro);
      }
      else {
        std::cout << " Extension to macro element not implemented for those element" << std::endl;
      }
    }



    // -----------------------------------------------------------------------------
    template<typename Mesh>
    void boundPreservingLimiter_P1(const FunFEM<Mesh>& uh, Rn& u_new, double min_val, double max_val, const MacroElement<Mesh>& macro) {
      typedef typename Mesh::Element Element;
      typedef typename Mesh::Rd Rd;
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename FElement::QF QF;
      const QF& qf(*QF_Simplex<typename FElement::RdHat>(2));

      u_new = uh.v;

      int nv_loc = Mesh::Rd::d +1;
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const ActiveMesh<Mesh>& Kh(macro.Th_);
      std::map<int, double>  map_mean_value;
      Rn uM(u_new);
      FunFEM<Mesh> fun_uM(Wh, uM);

      extendToMacro(uh, uM, map_mean_value, macro);
      u_new = uM;

      // Apply limiter in all element
      for(int k=Wh.first_element(); k<Wh.last_element();++k) {

        const FElement& FK(Wh[k]);
        const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
        double u_bar_K = 0.;
        double max_K=-1e300, min_K = 1e300;
        double areaCut = 0.;
        for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
          double meas = cutK.measure(it);
          areaCut += meas;
          // get the mean value on K
          for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip = cutK.mapToPhysicalElement(it, ip);
            double Cint = meas * ip.getWeight();
            // u_bar_K += Cint*uh.eval(k, mip);
            u_bar_K += Cint*fun_uM.eval(k, mip);
          }

          // get the local min and max on K
          for(int ipq = 0; ipq < nv_loc; ++ipq)  {
            Rd mip = cutK.get_vertex(it, ipq);
            // double val = uh.eval(k, mip);
            double val = fun_uM.eval(k, mip);

            max_K  = max(max_K, val);
            min_K  = min(min_K, val);
          }
        }
        u_bar_K = u_bar_K / areaCut;

        // 3) compute theta
        double v1 = fabs((min_val - u_bar_K) / (min_K - u_bar_K));
        if((min_K - u_bar_K) == 0.) v1 = 1.;//
        double v2 = fabs((max_val - u_bar_K) / (max_K - u_bar_K));
        if((max_K - u_bar_K) == 0.) v2 = 1.;
        double theta = min(min(v1, v2), 1.);
        // if( fabs(theta) < Epsilon) theta = 0.;

        // 4) replace the dof
        for(int df=0; df<FK.NbDoF() ; ++df) {
          int df_glb = FK.loc2glb(df);
          u_new(df_glb) = theta*(uh(df_glb)-u_bar_K) + u_bar_K;


          // if(u_new(df_glb) < min_val && u_new(df_glb) < -1e-6) {
          //   std::cout << " ----------------------- " << std::endl;
          //   std::cout << " Element " << k << "\n"
          //   << " theta   " << theta << "\n"
          //   << " min K   " << min_K << "\n"
          //   << " minK-u  " << min_K - u_bar_K << "\n"
          //   << " minv-u  " << min_val - u_bar_K << "\n"
          //   << " u_bar_K " << u_bar_K << "\n"
          //   << " u_M     " << uM(df_glb) << "\n"
          //   << " u_new   " << u_new(df_glb) << "\n"
          //   << " uh      " << uh(df_glb) << std::endl;
          // }


          // if(u_new(df_glb) - min_val < 0 && u_new(df_glb) - min_val > -Epsilon ) u_new(df_glb) = min_val;
          // if(u_new(df_glb) - max_val > 0 && u_new(df_glb) - max_val < Epsilon ) u_new(df_glb) = max_val;
          // if(fabs(u_new(df_glb) - min_val) < Epsilon ) u_new(df_glb) = min_val;
          // if(fabs(u_new(df_glb) - max_val) < Epsilon ) u_new(df_glb) = max_val;
        }
      }

      // loop over macro element
      for (auto& p : macro.macro_element) {
        const MElement& MK(p.second);
        int idx_root = MK.idx_root_element;

        // 1) the mean value
        double u_bar_M = map_mean_value[idx_root];
        // 2) get local min and max
        double max_M=-1e300, min_M = 1e300;
        for(auto k : MK.idx_element ) {
          const FElement& FK(Wh[k]);
          const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
          for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
            for(int ipq = 0; ipq < nv_loc; ++ipq)  {
              Rd mip = cutK.get_vertex(it, ipq);
              double val = fun_uM.eval(k, mip);
              max_M  = max(max_M, val);
              min_M  = min(min_M, val);
            }
          }
        }

        // 3) compute theta
        double v1 = fabs((min_val - u_bar_M) / (min_M - u_bar_M-Epsilon));
        double v2 = fabs((max_val - u_bar_M) / (max_M - u_bar_M+Epsilon));
        double theta = min(min(v1, v2), 1.);
        // std::cout << theta << std::endl;
        if( fabs(theta) < Epsilon) theta = 0.;
        // 4) replace the dof
        for(auto k : MK.idx_element ) {
          const FElement& FK(Wh[k]);
          for(int df=0; df<FK.NbDoF() ; ++df) {
            int df_glb = FK.loc2glb(df);
            u_new(df_glb) = theta*(uM(df_glb)-u_bar_M) + u_bar_M;

            // if(idx_root == 15545) {

            // if(u_new(df_glb) - min_val < 0 && u_new(df_glb) - min_val > -100*Epsilon ) u_new(df_glb) = min_val;
            // if(u_new(df_glb) - max_val > 0 && u_new(df_glb) - max_val < Epsilon )      u_new(df_glb) = max_val;
            // if(fabs(u_new(df_glb) - min_val) < Epsilon ) u_new(df_glb) = min_val;
            // if(fabs(u_new(df_glb) - max_val) < Epsilon ) u_new(df_glb) = max_val;

            // if(u_new(df_glb) < min_val && u_new(df_glb) < -1e-6) {
            //   std::cout << " ----------------------- " << std::endl;
            //   std::cout << " Element " << k << "\n"
            //             << " theta   " << theta << "\n"
            //             << " min K   " << min_M << "\n"
            //             << " minK-u  " << min_M - u_bar_M << "\n"
            //             << " minv-u  " << min_val - u_bar_M << "\n"
            //             << " u_bar_K " << u_bar_M << "\n"
            //             << " uh      " << uh(df_glb) << "\n"
            //             << " u_M     " << uM(df_glb) << "\n"
            //             << " u_new   " << u_new(df_glb) << std::endl;
            //             // getchar();
            // }
          }
        }
      }

    }


    template<typename Mesh>
    void applyBoundPreservingLimiter(const FunFEM<Mesh>& uh, Rn& u_new, double min_val, double max_val, const MacroElement<Mesh>& macro) {
      if(uh.getBasisFctType() == BasisFctType::P1dc) {
        boundPreservingLimiter_P1(uh, u_new, min_val, max_val, macro);
      }
      else {
        std::cout << " No boundary preserving limiter implemented for this elements" << std::endl;
      }
    }


    template<typename Mesh>
    void check_mean_value(const FunFEM<Mesh>& uh, const MacroElement<Mesh>& macro, double minV, double maxV, Rn& u_mean) {
      typedef typename Mesh::Element Element;
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename FElement::QF QF;
      typedef typename Mesh::Rd Rd;

      const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
      const GFESpace<Mesh>& Wh(*uh.Vh);

      if(u_mean.size() != Wh.get_nb_element()){
        u_mean.resize(Wh.get_nb_element());
      }
      u_mean = 0.;

      // compute  on each element
      for(int k=Wh.first_element(); k<Wh.last_element();++k) {

        const FElement& FK(Wh[k]);
        const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
        double u_bar_K = 0.;
        double areaCut = 0.;
        for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
          double meas = cutK.measure(it);
          areaCut += meas;
          // get the mean value on K
          for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip = cutK.mapToPhysicalElement(it, ip);
            double Cint = meas * ip.getWeight();
            u_bar_K += Cint*uh.eval(k, mip);
          }
        }
        u_bar_K = u_bar_K / areaCut;
        u_mean[k] = u_bar_K;

      }

      // compute mean value of the macro elements
      for (auto& p : macro.macro_element) {
        const MElement& MK(p.second);
        int idx_root = MK.idx_root_element;
        double area_total = MK.area_total_;
        double areaCut = 0.;
        double u_bar_K = 0.;
        for(auto k : MK.idx_element ) {
          const FElement& FK(Wh[k]);
          const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
          for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
            double meas = cutK.measure(it);
            areaCut += meas;
            for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
              typename QF::QuadraturePoint ip(qf[ipq]); // integration point
              Rd mip = cutK.mapToPhysicalElement(it, ip);
              double Cint = meas * ip.getWeight();
              u_bar_K += Cint*uh.eval(k, mip);
            }
          }
        }
        assert(fabs(area_total - areaCut) < 1e-14);
        u_bar_K = u_bar_K / areaCut;
        for(auto k : MK.idx_element ) {
          u_mean[k] = u_bar_K;
        }
      }
      for(int i=0;i<u_mean.size();++i){
        double u_bar_K = u_mean(i);
        if(u_bar_K < minV-Epsilon || maxV+Epsilon < u_bar_K){
          std::cout << " Mean value does not satisfies maximum principle " <<  std::endl;
          std::cout << "[m, M] = [ " << minV << "\t" << u_bar_K << "\t"  << maxV << " ]" << std::endl;
        }
      }
    }
    // template<typename Mesh>
    // void check_mean_value(const FunFEM<Mesh>& uh, const MacroElement<Mesh>& macro, double minV, double maxV, Rn& u_mean) {
    //   typedef typename Mesh::Element Element;
    //   typedef typename GFESpace<Mesh>::FElement FElement;
    //   typedef typename FElement::QF QF;
    //   typedef typename Mesh::Rd Rd;
    //
    //   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
    //   const GFESpace<Mesh>& Wh(*uh.Vh);
    //
    //   if(u_mean.size() != Wh.get_nb_element()){
    //     u_mean.resize(Wh.get_nb_element());
    //   }
    //   u_mean = 0.;
    //
    //   // compute mean value on each element
    //   for(int k=Wh.first_element(); k<Wh.last_element();++k) {
    //
    //     const FElement& FK(Wh[k]);
    //     const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
    //     double u_bar_K = 0.;
    //     double areaCut = 0.;
    //     for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
    //       double meas = cutK.measure(it);
    //       areaCut += meas;
    //       // get the mean value on K
    //       for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
    //         typename QF::QuadraturePoint ip(qf[ipq]); // integration point
    //         Rd mip = cutK.mapToPhysicalElement(it, ip);
    //         double Cint = meas * ip.getWeight();
    //         u_bar_K += Cint*uh.eval(k, mip);
    //       }
    //     }
    //     u_bar_K = u_bar_K / areaCut;
    //     u_mean[k] = u_bar_K;
    //   }
    //
    //   // compute mean value of the macro elements
    //   for (auto& p : macro.macro_element) {
    //     const MElement& MK(p.second);
    //     int idx_root = MK.idx_root_element;
    //     double area_total = MK.area_total_;
    //     double areaCut = 0.;
    //     double u_bar_K = 0.;
    //     for(auto k : MK.idx_element ) {
    //       const FElement& FK(Wh[k]);
    //       const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
    //       for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
    //         double meas = cutK.measure(it);
    //         areaCut += meas;
    //         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
    //           typename QF::QuadraturePoint ip(qf[ipq]); // integration point
    //           Rd mip = cutK.mapToPhysicalElement(it, ip);
    //           double Cint = meas * ip.getWeight();
    //           u_bar_K += Cint*uh.eval(k, mip);
    //         }
    //       }
    //     }
    //     assert(fabs(area_total - areaCut) < 1e-14);
    //     u_bar_K = u_bar_K / areaCut;
    //   }
    // }


  }

  namespace FEM {

    template<typename Mesh>
    void check_mean_value(const FunFEM<Mesh>& uh, double minV, double maxV, Rn& u_mean) {
      typedef typename Mesh::Element Element;
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename FElement::QF QF;
      typedef typename Mesh::Rd Rd;

      const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const Mesh& Kh(Wh.Th);

      if(u_mean.size() != Wh.get_nb_element()){
        u_mean.resize(Wh.get_nb_element());
      }
      u_mean = 0.;


      for(int k=Wh.first_element(); k<Wh.last_element();++k) {

        const FElement& FK(Wh[k]);
        const Element& K(FK.T);

        double u_bar_K = 0.;
        double meas = K.measure();
        // get the mean value on K
        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
          typename QF::QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = K.mapToPhysicalElement(ip);
          double Cint = ip.getWeight();
          u_bar_K += Cint*uh.eval(k, mip);
        }
        u_mean[k] = u_bar_K;
        if(u_bar_K < minV || maxV < u_bar_K){
          std::cout << " Mean value does not satisfies maximum principle " <<  std::endl;
          std::cout << "[m, M] = [ " << minV << "\t" << u_bar_K << "\t"  << maxV << " ]" << std::endl;

        }
      }
    }

    template<typename Mesh>
    void check_mean_value(const FunFEM<Mesh>& uh, double minV, double maxV) {
      typedef typename Mesh::Element Element;
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename FElement::QF QF;
      typedef typename Mesh::Rd Rd;

      const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const Mesh& Kh(Wh.Th);


      for(int k=Wh.first_element(); k<Wh.last_element();++k) {

        const FElement& FK(Wh[k]);
        const Element& K(FK.T);

        double u_bar_K = 0.;
        double meas = K.measure();
        // get the mean value on K
        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
          typename QF::QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = K.mapToPhysicalElement(ip);
          double Cint = ip.getWeight();
          u_bar_K += Cint*uh.eval(k, mip);
        }
        if(u_bar_K < minV || maxV < u_bar_K){
          std::cout << " Mean value does not satisfies maximum principle " <<  std::endl;
          std::cout << "[m, M] = [ " << minV << "\t" << u_bar_K << "\t"  << maxV << " ]" << std::endl;
          getchar();
        }
      }
    }


    template<typename Mesh>
    void limiter_Pei_P1(const FunFEM<Mesh>& uh, Rn& u_new, double min_val, double max_val) {
      typedef typename Mesh::Element Element;
      typedef typename Mesh::Rd Rd;
      typedef typename GFESpace<Mesh>::FElement FElement;
      typedef typename FElement::QF QF;
      const QF& qf(*QF_Simplex<typename FElement::RdHat>(2));

      u_new = uh.v;
      int nv_loc = Mesh::Rd::d +1;
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const Mesh& Kh(Wh.Th);

      // Apply limiter in all element
      for(int k=Wh.first_element(); k<Wh.last_element();++k) {

        const FElement& FK(Wh[k]);
        const Element& K(FK.T);

        double u_bar_K = 0.;
        double max_K=-1e300, min_K = 1e300;
        double meas = K.measure();
        // get the mean value on K
        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
          typename QF::QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = K.mapToPhysicalElement(ip);
          double Cint = ip.getWeight();
          u_bar_K += Cint*uh.eval(k, mip);
        }
        // if(u_bar_K < min_val ) u_bar_K = min_val;

        // get the local min and max on K
        for(int ipq = 0; ipq < nv_loc; ++ipq)  {
          Rd mip = K[ipq];
          double val = uh.eval(k, mip);

          max_K  = max(max_K, val);
          min_K  = min(min_K, val);
        }


        // 3) compute theta
        double v1 = fabs((min_val - u_bar_K) / (min_K - u_bar_K));
        if((min_K - u_bar_K) == 0.) v1 = 1.;//
        double v2 = fabs((max_val - u_bar_K) / (max_K - u_bar_K));
        if((max_K - u_bar_K) == 0.) v2 = 1.;
        double theta = min(min(v1, v2), 1.);

        // if( fabs(theta) < Epsilon) theta = 0.;

        // 4) replace the dof
        for(int df=0; df<FK.NbDoF() ; ++df) {
          int df_glb = FK.loc2glb(df);
          u_new(df_glb) = theta*(uh(df_glb)-u_bar_K) + u_bar_K;

          // if(u_new(df_glb) - min_val < 0 && u_new(df_glb) - min_val > -Epsilon ) u_new(df_glb) = min_val;
          // if(u_new(df_glb) - max_val > 0 && u_new(df_glb) - max_val < Epsilon ) u_new(df_glb) = max_val;

          // if(u_new(df_glb) < min_val) {
          //   std::cout << " ----------------------- " << std::endl;
          //   std::cout << " Element " << k << "\n"
          //             << " theta   " << theta << "\n"
          //             << " min K   " << min_K << "\n"
          //             << " minK-u  " << min_K - u_bar_K << "\n"
          //             << " minv-u  " << min_val - u_bar_K << "\n"
          //             << " u_bar_K " << u_bar_K << "\n"
          //             << " uh      " << uh(df_glb) << "\n"
          //             << " u_new   " << u_new(df_glb) << std::endl;
          //             getchar();
          //
          // }

        }
      }

    }

    template<typename Mesh>
    void minmaxP1(const FunFEM<Mesh>& uh, double& min_val, double& max_val) {
      typedef typename Mesh::Rd Rd;
      typedef typename Mesh::Element Element;
      typedef typename GFESpace<Mesh>::FElement FElement;
      int nv_loc = Mesh::Rd::d +1;
      const GFESpace<Mesh>& Wh(*uh.Vh);
      const Mesh& Kh(Wh.Th);

      min_val=1e300;
      max_val = -1e300;
      // Apply limiter in all element
      for(int k=Wh.first_element(); k<Wh.last_element();++k) {
        const FElement& FK(Wh[k]);

        // get the local min and max on K
        for(int ipq = 0; ipq < nv_loc; ++ipq)  {
          Rd mip = FK.T[ipq];
          // double val = uh.eval(k, mip);
          double val = uh.eval(k, mip);
          max_val  = max(max_val, val);
          min_val  = min(min_val, val);
        }
      }
    }

  };


};

// class Limiter {
//
//   typedef typename Mesh2::Element Element;
//   typedef typename FESpace2::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::QFB QFB;
//
//   const QF&  qf ;
//   const QFB& qfb;
//
// public:
//   Rn indicator;
//
//   Limiter() : qf(*QF_Simplex<R2>(3)), qfb(*QF_Simplex<R1>(3)){}
//
//   void KXRCF_indicator(const Fun2_h& uh, const Fun2_h& flux);
//
//   void limiting(const Fun2_h& uh, Rn& lh);
//   void min_and_max_average_neighbor(const FElement& FK, const Fun2_h&uh, double&m, double&M);
//   double average_on_Kj(const FElement& FK, const Fun2_h&uh);
//   double compute_alpha(const FElement& FK, const Fun2_h&uh, double Uj, double m, double M);
// };
//
//
// class Box_Side {
// public:
//   // a side is the part that satisfies ax+by+c < 0
//   double a, b, c;
//   Box_Side(double aa, double bb, double cc) : a(aa), b(bb), c(cc) {}
//
//   double f(const R2 A) const {
//     return (a*A.x+b*A.y+c);
//   }
//   bool is_negative(const R2 P) const {
//     return (f(P) <=0);
//   }
//   bool is_positive(const R2 P) const {
//     return (f(P) > 0);
//   }
//   bool cut_edge(const R2 A, const R2 B) const {
//     return (f(A)*f(B)< 0);
//   }
//   R2 get_cut_node(const R2 A, const R2 B) const {
//     double t = -f(A)/(f(B)-f(A));
//     return (1-t) * A + t * B;
//   }
//
// };
//
// class Filter_Box {
//   typedef typename Mesh2::Element Element;
// public:
//   static const int nv = 3;                         // triangle
//   typedef SortArray<Ubyte, nv> TriaIdx;           // the vertices of a triangle of the cut:
//
//   std::vector<Box_Side> sides_;
//   std::vector<R2>      vertices_;
//   std::vector<TriaIdx> triangles_;               // idx of the vertices
//
//   R2 minXY, maxXY;
//
//   Filter_Box(const R2 mminXY, const R2 mmaxXY) : minXY(mminXY), maxXY(mmaxXY){
//     sides_.push_back(Box_Side( 0., -1.,  minXY.y));
//     sides_.push_back(Box_Side( 1.,  0., -maxXY.x));
//     sides_.push_back(Box_Side( 0.,  1., -maxXY.y));
//     sides_.push_back(Box_Side(-1.,  0.,  minXY.x));
//   }
//
//   bool contain(const R2 P) const{
//     for(auto it=sides_.begin(); it != sides_.end();++it) {
//       if(it->is_positive(P)) return false;
//     }
//     return true;
//   }
//
//   bool intersect_with(const typename Mesh2::Element& K) const{
//     // if we find one node inside the box it is enough to know
//     for(int i=0; i<3;++i) {
//       if(this->contain((R2) K[i])) return true;
//     }
//     return false;
//   }
//
//   void cut_triangle(const typename Mesh2::Element& K) ;
//
//   void print(std::string filename = "newTriangle.dat") const ;
//
// };
//
// class Filter {
// public:
//
//
// };


#endif
