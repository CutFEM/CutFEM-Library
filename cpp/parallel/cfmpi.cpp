/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#include "cfmpi.hpp"

#if USE_MPI

//----------------- Static Members ---------------------------
int MPIcf::my_rank_                  = 0;
int MPIcf::size_                     = 1;
bool MPIcf::usePetsc_                = true;
MuteStdOstream *MPIcf::mute_         = 0;
MPI_Request *MPIcf::rq               = 0;
int MPIcf::loopDivision              = 0;
const MPI_Comm &MPIcf::Communicator_ = MPI_COMM_WORLD;

// constructor of the class
MPIcf::MPIcf(int &argc, char **&argv) {
    MPI_Init(&argc, &argv);                   // initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_); // set the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size_);    // set the size

    mute_ = new MuteStdOstream(); // mute ostream for non master
    muteStdOstreams();

    std::cout << "init parallel rank " << my_rank_ << " of " << size_ << std::endl;
    // PetscInitialize(&argc, &argv, NULL, NULL);
}

MPIcf::MPIcf() {

    MPI_Init(nullptr, nullptr); // initialize MPI

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_); // set the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size_);    // set the size

    mute_ = new MuteStdOstream(); // mute ostream for non master
    muteStdOstreams();
    std::cout << "init parallel rank " << my_rank_ << " of " << size_ << std::endl;

    usePetsc_ = false;
}

MPIcf::~MPIcf() {

    // if(usePetsc_)PetscFinalize();

    MPI_Finalize();
    std::cout << " \n MPI finalized correctly \n" << std::flush;

    size_ = 0;
    RecoverStdOstreams();
    delete mute_;
}

#endif