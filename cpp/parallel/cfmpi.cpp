#include "cfmpi.hpp"
// #include "cutFEMConfig.h"

//#include "petsc.h"

//----------------- Static Members ---------------------------
int MPIcf::my_rank_ = 0;
int MPIcf::size_ = 1;
bool MPIcf::usePetsc_ = true;
MuteStdOstream*  MPIcf::mute_ = 0;
MPI_Request* MPIcf::rq = 0;
int MPIcf::loopDivision = 0;
const MPI_Comm& MPIcf::Communicator_ = MPI_COMM_WORLD;

// constructor of the class
MPIcf::MPIcf(int &argc, char **& argv) {
  MPI_Init(&argc,&argv);                              // initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);           // set the rank
  MPI_Comm_size(MPI_COMM_WORLD, &size_);              // set the size

  mute_ = new MuteStdOstream();                       // mute ostream for non master
  muteStdOstreams();

  std::cout << "init parallele rank " <<  my_rank_ << " on " << size_ << std::endl;
  //PetscInitialize(&argc, &argv, NULL, NULL);

}

MPIcf::MPIcf() {

  MPI_Init(nullptr, nullptr);                                         // initialize MPI

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);           // set the rank
  MPI_Comm_size(MPI_COMM_WORLD, &size_);              // set the size

  mute_ = new MuteStdOstream();                       // mute ostream for non master
  muteStdOstreams();
  std::cout << "init parallele rank " <<  my_rank_ << " on " << size_ << std::endl;

  usePetsc_ = false;
}

MPIcf::~MPIcf() {

  // if(usePetsc_)PetscFinalize();

  MPI_Finalize();
  std::cout << " \n MPI finalize correctly \n" << std::flush ;

  size_=0;
  RecoverStdOstreams();
  delete mute_;

}
