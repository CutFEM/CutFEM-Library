#ifndef TIME_STUFF_HPP_
#define TIME_STUFF_HPP_

template<typename Mesh>
void set_velocity(FunFEM<Mesh>& uh, FunFEM<Mesh>& vel, FunFEM<Mesh>& ls, R tt)  {

  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Partition Partition;
  typedef typename Mesh::Rd Rd;

  const FESpace& velVh(*vel.Vh);
  const FESpace& Vh(*uh.Vh);
  vel.v = 0.;
  double loc_ls[Mesh::Element::nv];

  for(int k=0;k<velVh.NbElement();++k) {
    const FElement& FKvel((velVh)[k]);
    int kback = velVh.Th(FKvel.T);

    ls.eval(loc_ls, kback);
    const Partition& cutK =  Partition(FKvel.T, loc_ls);

    if(cutK.is_cut()) {
      int kcut0 = Vh.idxElementFromBackMesh(kback,0);
      int kcut1 = Vh.idxElementFromBackMesh(kback,1);

      for(int i=FKvel.dfcbegin(0);i<FKvel.dfcend(0);++i) {
        Rd x = FKvel.Pt(i);
        R lsval = ls.eval(kback, x);
        int kcut = (lsval > 0)? kcut0 : kcut1;
        for(int ci=0;ci<Rd::d;++ci) {
          vel( FKvel(i + FKvel.dfcbegin(ci))) = uh.eval(kcut, x, tt, ci, op_id);
        }
      }
    }
    else{
      int k = (loc_ls[0] > 0)? Vh.idxElementFromBackMesh(kback,0)
      : Vh.idxElementFromBackMesh(kback,1);

      for(int i=FKvel.dfcbegin(0);i<FKvel.dfcend(0);++i) {
        Rd x = FKvel.Pt(i);
        for(int ci=0;ci<Rd::d;++ci) {
          vel( FKvel(i + FKvel.dfcbegin(ci))) = uh.eval(k, x, tt, ci, op_id);
        }
      }
    }
  }
}




#endif
