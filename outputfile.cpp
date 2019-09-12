// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include "outputfile.hpp"
#include "parallelenvironment.hpp"

output_file::output_file(const std::string& file_name, long chunk_size){
	// Bug in MPI
  MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(file_name.c_str()), MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &f);
  MPI_File_close(&f);
  //
  MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(file_name.c_str()), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &f);
  long chunk_sizes_in[parallel_environment::get_num_processes()];
  long chunk_sizes_out[parallel_environment::get_num_processes()];
  int i;
  for(i=0; i<parallel_environment::get_num_processes(); ++i){
	  chunk_sizes_out[i]=chunk_size;
  }
  MPI_Alltoall(chunk_sizes_out, 1, MPI_LONG, chunk_sizes_in, 1, MPI_LONG, MPI_COMM_WORLD);
 
  offset=0;
  for(i=0; i<parallel_environment::get_process_num(); ++i){
    offset=offset+chunk_sizes_in[i];
  }
  MPI_File_set_view(f, offset, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, "native", MPI_INFO_NULL);
}

output_file::~output_file(){
	MPI_File_close(&f);
}

int output_file::write(const std::string& text){
   MPI_Status status;
   return MPI_File_write(f, const_cast<char*>(text.c_str()), text.length(), MPI_UNSIGNED_CHAR, &status);
}
