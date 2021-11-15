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
