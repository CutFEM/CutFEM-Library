
template<typename M>
FunFEM<M>::FunFEM(const FESpace& vh, const ExpressionVirtual& fh) :
FunFEMVirtual (vh.NbDoF()),
alloc(true), Vh(&vh),
// data(new double[vh.NbDoF()]),
// v(data, vh.NbDoF()) ,
databf(new double[10*vh[0].NbDoF()*vh.N*4])
{
  assert(Vh->N == 1);
  double dataSend[Vh->nbDoF];
  Rn_ fhSend(dataSend,Vh->nbDoF); fhSend = 1e+50;
  const int d = Vh->N;
  const int nve = Vh->TFE(0)->NbPtforInterpolation;
  KNM<R>   Vpf(nve,d);                       // value of f at the interpolation points
  KN<R> ggf(Vh->MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    const FElement& FK((*Vh)[k]);
    const int nbdf = FK.NbDoF();            // nof local
    const int domain = FK.whichDomain();
    const int kb = Vh->idxElementInBackMesh(k);

    for (int p=0;p<FK.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(FK.Pt(p));       // the coordinate of P in K hat
      for(int i=0;i<d;++i) {
        Vpf(p,i) = fh.evalOnBackMesh(kb, domain, P);
      }
    }
    std::cout << Vpf << std::endl;
    FK.Pi_h(Vpf,ggf);
    for (int df=0;df<nbdf;df++) {
      fhSend(FK.loc2glb(df)) = ggf[df] ;
      // fh[K(df)] =  ggf[df] ;
    }
    // for(int j=FK.dfcbegin(0);j<FK.dfcend(0);++j) {
    //   Rd mip = FK.Pt(j);
    //   fhSend(FK.loc2glb(j)) = fh.evalOnBackMesh(kb, domain, mip);
    //   // v(FK.loc2glb(j)) = fh.evalOnBackMesh(kb, domain, mip);
    // }
    getchar();
  }
  MPIcf::AllReduce(dataSend, data, fhSend.size(),MPI_MIN);
}

template<typename M>
FunFEM<M>::FunFEM(const FESpace& vh, const ExpressionVirtual& fh1, const ExpressionVirtual& fh2) :
FunFEMVirtual (vh.NbDoF()),
alloc(true), Vh(&vh),
// data(new double[vh.NbDoF()]),
// v(data, vh.NbDoF()) ,
databf(new double[10*vh[0].NbDoF()*vh.N*4])
{
  assert(Vh->N == 2);
  double dataSend[Vh->nbDoF];
  Rn_ fhSend(dataSend,Vh->nbDoF); fhSend = 1e+50;
  const int d = Vh->N;
  const int nve = Vh->TFE(0)->NbPtforInterpolation;
  KNM<R>   Vpf(nve,d);                       // value of f at the interpolation points
  KN<R> ggf(Vh->MaxNbDFPerElement);           // stock the values of the dof of the interpolate

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    const FElement& FK((*Vh)[k]);
    const int nbdf = FK.NbDoF();            // nof local
    const int domain = FK.whichDomain();
    const int kb = Vh->idxElementInBackMesh(k);


    for (int p=0;p<FK.tfe->NbPtforInterpolation;p++) {      // all interpolation points
      const Rd & P(FK.Pt(p));       // the coordinate of P in K hat
      for(int i=0;i<d;++i) {
        const ExpressionVirtual& fh = (d==0)?fh1 : fh2;
        Vpf(p,i) = fh.evalOnBackMesh(kb, domain, P);
      }
    }
    // std::cout << Vpf << std::endl;
    FK.Pi_h(Vpf,ggf);
    for (int df=0;df<nbdf;df++) {
      fhSend(FK.loc2glb(df)) = ggf[df] ;
      // fh[K(df)] =  ggf[df] ;
    }
    // for(int i=0, j=FK.dfcbegin(ci);j<FK.dfcend(ci);++j,++i) {
    //   Rd mip = FK.Pt(i);
    //   fhSend(FK.loc2glb(j)) = fh.evalOnBackMesh(kb, domain, mip);
    // }

  }
  MPIcf::AllReduce(dataSend, data, fhSend.size(),MPI_MIN);
}



template<typename M>
void FunFEM<M>::print() const{
  std::cout << v << std::endl;
}

template<typename M>
double FunFEM<M>::eval(const int k, const R* x, int cu, int op) const{
  const FElement& FK((*Vh)[k]);
  int ndf = FK.NbDoF();
  RNMK_ w(databf, ndf,Vh->N,op_dz+1);
  FK.BF(FK.T.toKref(x), w);
  double val = 0.;

  for(int j=FK.dfcbegin(cu);j<FK.dfcend(cu);++j) {
     val += v[FK(j)]*w(j,cu,op);
  }
  return val;
}

template<typename M>
double FunFEM<M>::eval(const int k, const R* x, const R t, int cu, int op, int opt) const{

  if(!In) return eval(k,x,cu,op);

  const FElement& FK((*Vh)[k]);
  int ndf = FK.NbDoF();
  RNMK_ w(databf, ndf,Vh->N,op_dz+1);
  KNMK<R> wt(In->NbDoF(),1,op_dz);

  FK.BF(FK.T.toKref(x), w);
  In->BF(In->T.toKref(t), wt );


  double val = 0.;
  for(int jt=0; jt<In->NbDoF();++jt) {
    for(int j=FK.dfcbegin(cu);j<FK.dfcend(cu);++j){
      val += v[FK(j)+jt*Vh->NbDoF()]*w(j,cu,op)*wt(jt,0,opt);
    }
  }

  return val;
}

template<typename M>
void FunFEM<M>::eval(R* u, const int k) const {
  assert(v && u);
  const FElement& FK((*Vh)[k]);
  for(int ci=0; ci<Vh->N;++ci) {              // loop over componant
    for(int j=FK.dfcbegin(ci);j<FK.dfcend(ci);++j)
    u[j] = v[FK(j)];
  }
}

template<typename M>
double FunFEM<M>::evalOnBackMesh(const int kb, const R* x, int cu, int op, int dom) const{
  int k = Vh->idxElementFromBackMesh(kb,dom);
  return eval(k, x, cu, op);
}

template<typename M>
double FunFEM<M>::evalOnBackMesh(const int kb, const R* x, const R t, int cu, int op, int opt, int dom) const{

  int k = Vh->idxElementFromBackMesh(kb,dom);

  return eval(k, x, t, cu, op, opt);
}

template<typename M>
std::list<ExpressionFunFEM<M>> FunFEM<M>::expression(int n)const{
  if(n == -1) n = Vh->N;
  assert(n<= Vh->N);
  std::list<ExpressionFunFEM<Mesh>> l;
  for(int i=0;i<n;++i){
    // const ExpressionFunFEM<Mesh> e(*this, i, op_id);
    l.push_back(ExpressionFunFEM<Mesh>(*this, i, op_id));
  }
  return l;
}

// template<typename M>
// std::list<ExpressionFunFEM<M>> LaxFriedrichs_flux(const std::list<ExpressionFunFEM<M>>& U, int cu)
// {
//   std::list<ExpressionFunFEM<M>> tt;
//   auto it=U.begin();
//   std::list<ExpressionFunFEM<M>> l;
//   const ExpressionFunFEM<M>& rho(*it);
//   l.push_back(rho*rho);
//   l.push_back(rho*rho);
//   return l;
//   // auto it=U.begin();
//   // if(cu == 0) {
//   //   return {*(it+1) , *(it+2)};
//   // }
// }
