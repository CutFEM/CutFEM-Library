#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

#include "../common/RNM.hpp"
#include "FESpace.hpp"
//#include "../parallel/cfmpi.hpp"




/*
  Interpolate f : Rd->R    on space Vh
 - output : fh contains the values
 */
template<typename F>
void interpolate(const F& Mh, KN<double>& fh, R(*f)(const typename F::Rd ) ){
  // std::cout << " need to double check this interpolate function and add MPI" << std::endl;
  // assert(0);
  typedef typename F::Rd Rd;
  typedef typename F::Element::RdHat RdHat;
  // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
  fh.init(Mh.nbDoF);
  const int d = 1;
  const int nve = Mh.MaxNbNodePerElement;
  KNM<R>   Vpf(nve,1);                       // value of f at the interpolation points
  KN<R> ggf(Mh.MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  // for (int t=Mh.first_element();t<Mh.last_element();
  //      t+= Mh.next_element()) {      // loop over element
  for (int t=0;t<Mh.NbElement();  t+=1) {
    typename F::FElement K(Mh[t]);
    const int nbdf = K.NbDoF();            // nof local

    for (int p=0;p<K.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(K.Pt(p));       // the coordinate of P in K hat
      Vpf(p,0) = f(P);

    }

    K.Pi_h(Vpf,ggf);
    for (int df=0;df<nbdf;df++) {
      // fhSend[K(df)] =  ggf[df] ;
      fh[K(df)] =  ggf[df] ;

    }
  }

  // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}


/*
  Interpolate f : Rd->R    on space Vh
 - output : fh contains the values
 */
template<typename F>
void interpolate(const F& Mh, KN_<double>& fh, R(*f)(const typename F::Rd, int i ) ){
  // std::cout << " need to double check this interpolate function and add MPI" << std::endl;
  // assert(0);
  typedef typename F::Rd Rd;
  typedef typename F::Element::RdHat RdHat;
  // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
  assert(fh.size() == Mh.nbDoF);
  // fh.init(Mh.nbDoF);
  const int d = Mh.N;
  const int nve = Mh.TFE(0)->NbPtforInterpolation;
  KNM<R>   Vpf(nve,d);                       // value of f at the interpolation points
  KN<R> ggf(Mh.MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  // for (int t=Mh.first_element();t<Mh.last_element();
  //      t+= Mh.next_element()) {      // loop over element
  for (int t=0;t<Mh.NbElement();  t+=1) {
    typename F::FElement K(Mh[t]);
    const int nbdf = K.NbDoF();            // nof local

    for (int p=0;p<K.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(K.Pt(p));       // the coordinate of P in K hat
      for(int i=0;i<d;++i) {
        Vpf(p,i) = f(P,i);
      }
    }

    K.Pi_h(Vpf,ggf);
    for (int df=0;df<nbdf;df++) {
      // fhSend[K(df)] =  ggf[df] ;
      fh[K(df)] =  ggf[df] ;

    }
  }
  // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}


/*
  Interpolate f : Rd->R    on space Vh
 - output : fh contains the values
 */
template<typename F>
void interpolate(const F& Mh, KN_<double>& fh, R(*f)(const typename F::Rd, int ii, int dom ) ){
  // std::cout << " need to double check this interpolate function and add MPI" << std::endl;
  // assert(0);
  typedef typename F::Rd Rd;
  typedef typename F::Element::RdHat RdHat;
  // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
  assert(fh.size() == Mh.nbDoF);
  // fh.init(Mh.nbDoF);
  const int d = Mh.N;
  const int nve = Mh.TFE(0)->NbPtforInterpolation;
  KNM<R>   Vpf(nve,d);                       // value of f at the interpolation points
  KN<R> ggf(Mh.MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  // for (int t=Mh.first_element();t<Mh.last_element();
  //      t+= Mh.next_element()) {      // loop over element
  for (int t=0;t<Mh.NbElement();  t+=1) {
    typename F::FElement K(Mh[t]);
    const int nbdf = K.NbDoF();            // nof local
    const int domain = K.whichDomain();

    for (int p=0;p<K.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(K.Pt(p));       // the coordinate of P in K hat
      for(int i=0;i<d;++i) {
        Vpf(p,i) = f(P,i, domain);
      }
    }

    K.Pi_h(Vpf,ggf);

    for (int df=0;df<nbdf;df++) {
      fh[K(df)] =  ggf[df] ;
    }
  }

// MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}

/*
  Interpolate f : Rd->R    on space Vh
 - output : fh contains the values
 */
template<typename F>
void interpolate(const F& Mh, KN_<double>& fh, R(*f)(const typename F::Rd, int i, int domain, R t ), R tid ){
  // std::cout << " need to double check this interpolate function and add MPI" << std::endl;
  // assert(0);
  typedef typename F::Rd Rd;
  typedef typename F::Element::RdHat RdHat;
  // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
  assert(fh.size() == Mh.nbDoF);
  // fh.init(Mh.nbDoF);
  const int d = Mh.N;
  const int nve = Mh.TFE(0)->NbPtforInterpolation;
  KNM<R>   Vpf(nve,d);                       // value of f at the interpolation points
  KN<R> ggf(Mh.MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  // for (int t=Mh.first_element();t<Mh.last_element();
  //      t+= Mh.next_element()) {      // loop over element
  for (int t=0;t<Mh.NbElement();  t+=1) {
    typename F::FElement K(Mh[t]);
    const int nbdf = K.NbDoF();            // nof local
    const int domain = K.whichDomain();

    for (int p=0;p<K.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(K.Pt(p));       // the coordinate of P in K hat
      for(int i=0;i<d;++i) {
        Vpf(p,i) = f(P,i, domain, tid);
      }
    }

    K.Pi_h(Vpf,ggf);


    for (int df=0;df<nbdf;df++) {
      fh[K(df)] =  ggf[df] ;
    }
  }

  // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}



/*
  Interpolate f : Rd->R    on space Vh
 - output : fh contains the values
 */
template<typename F>
void interpolate(const F& Mh, KN_<double>& fh, R(*f)(const typename F::Rd, int, R ), R tid ){
  // std::cout << " need to double check this interpolate function and add MPI" << std::endl;
  // assert(0);
  typedef typename F::Rd Rd;
  typedef typename F::Element::RdHat RdHat;
  // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
  assert(fh.size() == Mh.nbDoF);
  // fh.init(Mh.nbDoF);
  const int d = Mh.N;
  const int nve = Mh.TFE(0)->NbPtforInterpolation;
  KNM<R>   Vpf(nve,d);                       // value of f at the interpolation points
  KN<R> ggf(Mh.MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  // for (int t=Mh.first_element();t<Mh.last_element();
  //      t+= Mh.next_element()) {      // loop over element
  for (int t=0;t<Mh.NbElement();  t+=1) {
    typename F::FElement K(Mh[t]);
    const int nbdf = K.NbDoF();            // nof local

    for (int p=0;p<K.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(K.Pt(p));       // the coordinate of P in K hat
      for(int i=0;i<d;++i) {
        Vpf(p,i) = f(P,i, tid);
      }
    }

    K.Pi_h(Vpf,ggf);

    for (int df=0;df<nbdf;df++) {
      // fhSend[K(df)] =  ggf[df] ;
      fh[K(df)] =  ggf[df] ;

    }
  }
  // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}



/*
  Interpolate f : Rd->R    on space time Vh
 - output : fh contains the values
 */
template<typename F>
void interpolate(const F& Mh, const TimeSlab& In, KN_<double>& fh, R(*f)(const typename F::Rd, int, R ) ){
  // std::cout << " need to double check this interpolate function and add MPI" << std::endl;
  // assert(0);
  typedef typename F::Rd Rd;
  typedef typename F::Element::RdHat RdHat;
  // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
  assert(fh.size() == Mh.nbDoF*In.NbDoF());
  // fh.init(Mh.nbDoF);
  const int d = Mh.N;
  const int nve = Mh.TFE(0)->NbPtforInterpolation;
  const int nvt = In.tfe->NbPtforInterpolation;
  KNM<R>   Vpt(nvt,1);
  KNM<R> ggft(Mh.MaxNbDFPerElement, In.NbDoF());           // stock the values of the dof of the interpolate


  KNMK<R>   Vpft(nve,d, nvt);                       // value of f at the interpolation points

  // for (int t=Mh.first_element();t<Mh.last_element();
  //      t+= Mh.next_element()) {      // loop over element
  for (int t=0;t<Mh.NbElement();  t+=1) {
    typename F::FElement K(Mh[t]);
    const int nbdf = K.NbDoF();            // nof local


    // compute all the value , space and time
    for(int it=0;it<In.tfe->NbPtforInterpolation;++it) {
      const R1 & tq(In.Pt(it));
      for (int p=0;p<K.tfe->NbPtforInterpolation;p++) {      // all interpolation points
        const Rd & P(K.Pt(p));       // the coordinate of P in K hat
        for(int i=0;i<d;++i) {
          Vpft(p,i,it) = f(P,i,tq);

        }
      }
    }

    // perfom the spacxe interpolation for each time dof
    for(int it=0;it<In.tfe->NbPtforInterpolation;++it) {
      KN_<R> ggg(ggft('.',it));
      K.Pi_h(Vpft('.','.',it),ggg);
    }

    // perfom time interpolation and save/replacve the value in the matrix by the good one
    for (int df=0;df<nbdf;df++) {
      for(int it=0;it<In.tfe->NbPtforInterpolation;++it) {
        Vpt(it,0) = ggft(df,it);
      }
      KN_<R> ggg(ggft(df,'.'));
      In.Pi_h(Vpt,ggg);
    }


    for(int it=0;it<In.tfe->NbPtforInterpolation;++it) {
      for (int df=0;df<nbdf;df++) {
        fh[K.loc2glb(df,it)] =  ggft(df,it);//[df] ;
      }
    }
  }
  // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}








#endif
