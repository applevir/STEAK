// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#ifndef OUTPUT_FILE
#define OUTPUT_FILE

#include <mpi.h>
#include <string>

class output_file{
private:
    int offset;
	MPI_File f;
public:
	output_file(const std::string& file_name, long chunk_size);
	~output_file();
	int write(const std::string& text);
};

#endif
