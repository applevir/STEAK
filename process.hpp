// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#ifndef PROCESS_HPP
#define PROCESS_HPP

#include "data.hpp"
#include <vector>
#include "ltr.hpp"
#include "smithwaterman.hpp"
#include <string>
#include "inputfile.hpp"
#include "buffer.hpp"

#define VECTOR_SIZE 2000

class match_result{
public:
   SmithWatermanResult swr;
   int trimmedPosition;
   std::string trimmedPhred, trimmedRead, matchRead, matchPhred;
   bool selected;
   match_result(){
	 selected=false;
   }
};

class result_texts{
public:
    std::string fastq1_txt, fastq2_txt, fastq_te_txt, read_info_txt;
};

int process_reads_paired(const double& alignment_quality,
		              const double& match_quality,  
					  int transposon_length, 
					  int host_length, 
					  const std::vector<LTR>& TE_references, 
                                          const std::vector<Read>& reads,
                                          std::string& fastq1,
                                          std::string& fastq2,
                                          std::string& fastq_te,
                                          std::string& read_info);


int process_reads_unpaired(const double& alignment_quality,
                              const double& match_quality,
                                          int transposon_length,
                                          int host_length,
                                          const std::vector<LTR>& TE_references,
                                          const std::vector<Read>& reads,
                                          std::string& fastq1,
                                          std::string& fastq_te,
                                          std::string& read_info);


int process_buffer_paired(const double& alignment_quality,
                              const double& match_quality,
                                          int transposon_length,
                                          int host_length,
                                          const std::vector<LTR>& TE_references,
			                  buffer<std::pair<Read, Read>>* b,
                                          result_texts* rt);



int process_buffer_unpaired(const double& alignment_quality,
                              const double& match_quality,
                                          int transposon_length,
                                          int host_length,
                                          const std::vector<LTR>& TE_references,
			                  buffer<Read>* b,
                                          std::string* fastq1_txt,
                                          std::string* fastq_te_txt,
                                          std::string* read_info_txt);



#endif
