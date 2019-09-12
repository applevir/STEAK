// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include "parallelenvironment.hpp"



class parallel_exception: public std::exception{
  virtual const char* what() const throw(){
    return "Unable to set up parallel environment";
  }
} parallel_exception;


parallel_environment::parallel_environment(int argc, char** argv){
  int rv;
  rv=MPI_Init(&argc, &argv);
  if(rv!=MPI_SUCCESS){
    MPI_Abort(MPI_COMM_WORLD, rv);
    throw parallel_exception;
  }
}

parallel_environment::~parallel_environment(){
  MPI_Finalize();
}
