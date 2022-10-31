#ifndef REINIT_HPP_
#define REINIT_HPP_

#include "baseCutProblem.hpp"
#include "projection.hpp"

template<class Mesh>
class CReinitialization {
  typedef GFESpace<Mesh> FESpace;
  typedef FunFEM<Mesh> Fun;
  typedef typename TypeInterface<Mesh::Rd::d>::Interface Interface;

  public :
  int number_iteration = 20;
  double epsilon_diffusion = 1e-3;
  string ON_OFF = "ON";
  string ODE_method = "Euler";
  double dt = 1e-2;
  string mass_correction = "ON";
  double precision_correction = 1e-5;
  int max_iteration_correction = 10;

  CReinitialization() {}

  void man() {
    std::cout << " Class Reinitialization \n";
    std::cout << " Fields of the class : " << std::endl;
    std::cout << " number_iteration         -> number of iteration of the reinitialization process \n"
              << " epsilon_diffusion        -> value of the diffusion epsilon of the streamLine diffusion method \n"
              << " ON_OFF                   -> allow to enable/disable the reinitialization (\"ON\" or \"OFF\") \n"
              << " ODE_method               -> method used to discretize the time derivative ( \"Euler\" or \"RK2\") \n"
              << " dt                       -> value of the time step use in the method \n"
              << " mass_correction          -> allow to enable/disable the correction of the mass (\"ON\" or \"OFF\") \n"
              << " precision_correction     -> precision of the mass correction  \n"
              << " max_iteration_correction -> maximun number of iteration if the precision is not reach \n"
              << std::endl;
            }
  void info() {
    std::cout << " Reinitialization : " << std::endl;
    std::cout << " number_iteration         -> " << number_iteration <<   "\n"
              << " epsilon_diffusion        -> " << epsilon_diffusion <<   "\n"
              << " ON_OFF                   -> " << ON_OFF <<   "\n"
              << " ODE_method               -> " << ODE_method  <<   "\n"
              << " dt                       -> " << dt <<   "\n"
              << " mass_correction          -> " << mass_correction << "\n"
              << " precision_correction     -> " << precision_correction <<   "\n"
              << " max_iteration_correction -> " << max_iteration_correction <<   "\n"
              << std::endl;
  }

  void perform(Fun& levelSet_P1);
  void perform(Fun& levelSet_k, Fun& levelSet_P1);

  double correction(Fun& levelSet, const double);
  double computeArea(const FESpace& Vh, const Interface& gamma);

};

template<class Mesh>
class Reinitialization : public ShapeOfLinProblem , public Solver{
  typedef GFESpace<Mesh> FESpace;
  typedef FunFEM<Mesh> Fun;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  typedef typename TypeInterface<Mesh::Rd::d>::Interface Interface;

  public :
  const FESpace& Vh;
  const Interface& interface;
  const QF& qf;


  const int number_iteration;
  const double epsilon_diffusion;
  const string ODE_method;
  const double dt;
  const string mass_correction;

  int ndofK;
  double* databf;

  Reinitialization(Fun& levelSet, const Interface& gam,
    int nnumber_iteration, double eepsilon_diffusion,
    string OODE_method, double ddt) :
  ShapeOfLinProblem(levelSet.Vh->NbDoF()),
  Solver(),
  Vh(*levelSet.Vh),
  interface(gam),
  qf(*QF_Simplex<Rd>(5)),
  number_iteration(nnumber_iteration),
  epsilon_diffusion(eepsilon_diffusion),
  ODE_method(OODE_method),
  dt(ddt)
  {
    ndofK = Vh[0].NbDoF(); // same local dof for every element
    databf = new double[Vh.NbDoF()*Vh.N*5];

    Fun u0(this->Vh);
    u0.init(levelSet.v);

    int iter = 0;

    while( iter < number_iteration ) {

      this->rhs = 0.0;
      this->assembly(levelSet, u0);

      Solver::solve(this->mat, this->rhs);

      Rn dw(this->rhs); dw -= levelSet.v;
      levelSet.init(this->rhs);
      iter ++;
    }


  }


  void assembly(const Fun& up, const Fun& u0);

};

#include "reinitialization.tpp"





#endif
