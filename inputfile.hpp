// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


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
      bool get_start_position(bool recursive=true);
      bool get_pair_beginning();
      bool get_read(Read& read, bool second_read=false);
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
