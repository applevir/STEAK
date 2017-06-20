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
