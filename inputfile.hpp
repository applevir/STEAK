#ifndef INPUTFILE_HPP
#define INPUTFILE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include "samdata.hpp"
#include "data.hpp"


class fasta_file{
  private:
    std::ifstream in;
    unsigned long long local_size, local_beginning, local_end, local_position;
    bool carry_on;
  public:
    fasta_file(const std::string filename, bool sequential=false);
    ~fasta_file();
    bool get_read(Read& read);
};


class fastq_file{
    private:
      std::ifstream in;
      unsigned long long local_size, local_beginning, local_end, local_position;
      bool carry_on;
    public:
      fastq_file(const std::string& filename);
      ~fastq_file();
      bool get_first_read(Read& read);
      bool get_read(Read& read);
};

class sam_file{
    private:
      std::ifstream in;
      unsigned long long local_size, local_beginning, local_end, local_position;
    public:
      sam_file(const std::string& filename);
      ~sam_file();
      bool get_read(Read& read);
};

class bam_pipe{
private:
	std::istream* input;
public:
	bam_pipe();
        int get_read(Read& read);
	
};


#endif
