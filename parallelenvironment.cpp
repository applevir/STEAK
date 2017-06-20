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
