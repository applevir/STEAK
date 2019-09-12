// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include "process.hpp"
#include "smithwaterman.hpp"
#include <map>
#include <boost/algorithm/string.hpp>
#include <omp.h>
using namespace std;



int process_reads_paired(const double& alignment_quality,
			 const double& match_quality,
			 int transposon_length,
			 int host_length,
			 const vector<LTR>& TE_references,
			 const vector<Read>& reads,
			 string& fastq1,
			 string& fastq2,
			 string& fastq_te,
			 string& read_info){

  if(reads.empty() || reads.size()>2){
    return -1;
  }

  int i;
  for(i=0; i<reads.size(); ++i){
    if(reads[i].cigarValue<=alignment_quality){
      break;
    }
  }

  if(i==reads.size()){
    return -1;
  }

  int i_ref, i_read;
  typedef pair<int, int> pair_t;
  map<pair<int, int>, match_result> results;
  int max_absolute_score=-1;
  pair<int, int> ref_read_pair_max;
  string matchRead;
  int best_TE_for_read[reads.size()], best_score_for_read[reads.size()];
  for(i=0; i<reads.size(); ++i){
    best_TE_for_read[i]=-1;
    best_score_for_read[i]=0;
  }

  for(i_ref=0; i_ref<TE_references.size(); ++i_ref){
    for(i_read=0; i_read<reads.size(); ++i_read){
#ifdef DEBUG
      std::cout << reads[i_read].name << "\n";
#endif
      smithwaterman(TE_references[i_ref].sequence, reads[i_read].sequence, results[pair_t(i_ref, i_read)].swr);
      if(results[pair_t(i_ref, i_read)].swr.relativeScore>=match_quality){
        getTrimmedValues(results[pair_t(i_ref, i_read)].swr, reads[i_read], TE_references[i_ref].direction, results[pair_t(i_ref, i_read)].trimmedRead, results[pair_t(i_ref, i_read)].trimmedPhred, results[pair_t(i_ref, i_read)].trimmedPosition, results[pair_t(i_ref, i_read)].matchRead, results[pair_t(i_ref, i_read)].matchPhred);
      }else{
	results[pair_t(i_ref, i_read)].trimmedRead=reads[i_read].sequence;
	results[pair_t(i_ref, i_read)].trimmedPhred=reads[i_read].phred;
	results[pair_t(i_ref, i_read)].trimmedPosition=0;
	results[pair_t(i_ref, i_read)].matchRead="";
	results[pair_t(i_ref, i_read)].matchPhred="";
      }
      if(results[pair_t(i_ref, i_read)].trimmedRead.length()>=transposon_length && results[pair_t(i_ref, i_read)].swr.relativeScore>=match_quality){
	results[pair_t(i_ref, i_read)].selected=true;
	if(results[pair_t(i_ref, i_read)].swr.score>max_absolute_score){
	  max_absolute_score=results[pair_t(i_ref, i_read)].swr.score;
	  ref_read_pair_max=pair_t(i_ref, i_read);
	}
	if(results[pair_t(i_ref, i_read)].swr.score>best_score_for_read[i_read]){
	  best_score_for_read[i_read]=results[pair_t(i_ref, i_read)].swr.score;
	  best_TE_for_read[i_read]=i_ref;
	}
      }
    }
  }

  if(max_absolute_score<0){
    return -1;
  }

  fastq1=       "@"+reads[ref_read_pair_max.second].name+"\n";
  fastq1=fastq1+results[ref_read_pair_max].trimmedRead+"\n";
  fastq1=fastq1+"+\n";
  fastq1=fastq1+results[ref_read_pair_max].trimmedPhred+"\n";
  int i_other_read=(ref_read_pair_max.second+1)%2;
  if(reads.size()>1){
    fastq2=       "@"+reads[i_other_read].name+"\n";
    if(best_TE_for_read[i_other_read]>=0){
      fastq2=fastq2+results[pair_t(best_TE_for_read[i_other_read], i_other_read)].trimmedRead+"\n";
      fastq2=fastq2+"+\n";
      fastq2=fastq2+results[pair_t(best_TE_for_read[i_other_read], i_other_read)].trimmedPhred+"\n";
    }else{
      fastq2=fastq2+reads[i_other_read].sequence+"\n";
      fastq2=fastq2+"+\n";
      fastq2=fastq2+reads[i_other_read].phred+"\n";
    }
  }

  fastq_te=         "@"+reads[ref_read_pair_max.second].name+"\n";
  fastq_te=fastq_te+results[ref_read_pair_max].matchRead+"\n";
  fastq_te=fastq_te+"+\n";
  fastq_te=fastq_te+results[ref_read_pair_max].matchPhred+"\n";

  string processedRead;
  matchRead=results[ref_read_pair_max].matchRead;
  boost::to_lower(matchRead);
  if(TE_references[ref_read_pair_max.first].direction<0){
    processedRead=results[ref_read_pair_max].trimmedRead+matchRead;
  }else{
    processedRead=matchRead+results[ref_read_pair_max].trimmedRead;
  }
  string processedMate;
  if(reads.size()>1){
    if(best_TE_for_read[i_other_read]>=0){
      matchRead=results[pair_t(best_TE_for_read[i_other_read], i_other_read)].matchRead;
      fastq_te=fastq_te+"@"+reads.front().name+"\n";
      fastq_te=fastq_te+matchRead+"\n";
      fastq_te=fastq_te+"+\n";
      fastq_te=fastq_te+results[pair_t(best_TE_for_read[i_other_read], i_other_read)].matchPhred+"\n";
      boost::to_lower(matchRead);
      if(TE_references[best_TE_for_read[i_other_read]].direction<0){
	processedMate=results[pair_t(best_TE_for_read[i_other_read], i_other_read)].trimmedRead+matchRead;
      }else{
	processedMate=matchRead+results[pair_t(best_TE_for_read[i_other_read], i_other_read)].trimmedRead;
      }
    }else{
      processedMate=reads[i_other_read].sequence;
    }
  }else{
    processedMate="NA";
  }
  read_info=reads.front().name+"\t"+TE_references.front().name+"\t"+to_string(results[ref_read_pair_max].trimmedRead.length())+"\t"+to_string(results[ref_read_pair_max].matchRead.length())+"\t"+to_string(results[ref_read_pair_max].swr.relativeScore)+"\t"+processedRead+"\t"+processedMate+"\t"+reads[ref_read_pair_max.second].chromosome+"|"+to_string(reads[ref_read_pair_max.second].position)+"|"+to_string(reads[ref_read_pair_max.second].quality)+"|"+reads[ref_read_pair_max.second].cigarString+"\t"+reads[i_other_read].chromosome+"|"+to_string(reads[i_other_read].position)+"|"+to_string(reads[i_other_read].quality)+"|"+reads[i_other_read].cigarString+"\n";
  return 0;
}



int process_reads_unpaired(const double& alignment_quality,
			   const double& match_quality,
			   int transposon_length,
			   int host_length,
			   const vector<LTR>& TE_references,
			   const vector<Read>& reads,
			   string& fastq1,
			   string& fastq_te,
			   string& read_info){


  if(reads.size()!=1){
    return -1;
  }

  int i;
  for(i=0; i<reads.size(); ++i){
    if(reads[i].cigarValue<=alignment_quality){
      break;
    }
  }

  if(i==reads.size()){
    return -1;
  }

  int i_ref, i_read;
  typedef pair<int, int> pair_t;
  map<pair<int, int>, match_result> results;
  int max_absolute_score=-1;
  pair<int, int> ref_read_pair_max;
  string matchRead;
  int best_TE_for_read[reads.size()], best_score_for_read[reads.size()];
  for(i=0; i<reads.size(); ++i){
    best_TE_for_read[i]=-1;
    best_score_for_read[i]=0;
  }

  for(i_ref=0; i_ref<TE_references.size(); ++i_ref){
    for(i_read=0; i_read<reads.size(); ++i_read){
#ifdef DEBUG
      std::cout << reads[i_read].name << "\n";
#endif
      smithwaterman(TE_references[i_ref].sequence, reads[i_read].sequence, results[pair_t(i_ref, i_read)].swr);
      if(results[pair_t(i_ref, i_read)].swr.relativeScore>=match_quality){
	getTrimmedValues(results[pair_t(i_ref, i_read)].swr, reads[i_read], TE_references[i_ref].direction, results[pair_t(i_ref, i_read)].trimmedRead, results[pair_t(i_ref, i_read)].trimmedPhred, results[pair_t(i_ref, i_read)].trimmedPosition, results[pair_t(i_ref, i_read)].matchRead, results[pair_t(i_ref, i_read)].matchPhred);
      }else{
	results[pair_t(i_ref, i_read)].trimmedRead=reads[i_read].sequence;
	results[pair_t(i_ref, i_read)].trimmedPhred=reads[i_read].phred;
	results[pair_t(i_ref, i_read)].trimmedPosition=0;
	results[pair_t(i_ref, i_read)].matchRead="";
	results[pair_t(i_ref, i_read)].matchPhred="";
      }
      if(results[pair_t(i_ref, i_read)].trimmedRead.length()>=transposon_length && results[pair_t(i_ref, i_read)].swr.relativeScore>=match_quality){
	results[pair_t(i_ref, i_read)].selected=true;
	if(results[pair_t(i_ref, i_read)].swr.score>max_absolute_score){
	  max_absolute_score=results[pair_t(i_ref, i_read)].swr.score;
	  ref_read_pair_max=pair_t(i_ref, i_read);
	}
	if(results[pair_t(i_ref, i_read)].swr.score>best_score_for_read[i_read]){
	  best_score_for_read[i_read]=results[pair_t(i_ref, i_read)].swr.score;
	  best_TE_for_read[i_read]=i_ref;
	}
      }
    }
  }

  if(max_absolute_score<0){
    return -1;
  }

  fastq1=       "@"+reads.front().name+"\n";
  fastq1=fastq1+results[ref_read_pair_max].trimmedRead+"\n";
  fastq1=fastq1+"+\n";
  fastq1=fastq1+results[ref_read_pair_max].trimmedPhred+"\n";

  fastq_te=         "@"+reads.front().name+"\n";
  fastq_te=fastq_te+results[ref_read_pair_max].matchRead+"\n";
  fastq_te=fastq_te+"+\n";
  fastq_te=fastq_te+results[ref_read_pair_max].matchPhred+"\n";

  string processedRead;
  matchRead=results[ref_read_pair_max].matchRead;
  boost::to_lower(matchRead);
  if(TE_references[ref_read_pair_max.first].direction<0){
    processedRead=results[ref_read_pair_max].trimmedRead+matchRead;
  }else{
    processedRead=matchRead+results[ref_read_pair_max].trimmedRead;
  }
  string processedMate;
  processedMate="NA";

  read_info=reads.front().name+"\t"+TE_references.front().name+"\t"+to_string(results[ref_read_pair_max].trimmedRead.length())+"\t"+to_string(results[ref_read_pair_max].matchRead.length())+"\t"+to_string(results[ref_read_pair_max].swr.relativeScore)+"\t"+processedRead+"\t"+processedMate+"\t"+reads[ref_read_pair_max.second].chromosome+"|"+to_string(reads[ref_read_pair_max.second].position)+"|"+to_string(reads[ref_read_pair_max.second].quality)+"|"+reads[ref_read_pair_max.second].cigarString+"\tNA\n";
  return 0;
}


int process_buffer_paired(const double &alignment_quality, const double &match_quality, int transposon_length, int host_length, const std::vector<LTR> &TE_references, buffer<std::pair<Read, Read>> *b, result_texts* rt){
  rt->fastq1_txt="";
  rt->fastq2_txt="";
  rt->fastq_te_txt="";
  rt->read_info_txt="";
  int i, j;
  int n_threads;
#pragma omp parallel
  n_threads=omp_get_num_threads();

  std::string fastq1_txt_array[n_threads];
  std::string fastq2_txt_array[n_threads];
  std::string fastq_te_txt_array[n_threads];
  std::string read_info_txt_array[n_threads];
  for(i=0; i<n_threads; ++i){
    fastq1_txt_array[i]="";
    fastq2_txt_array[i]="";
    fastq_te_txt_array[i]="";
    read_info_txt_array[i]="";
  }
  std::vector<std::pair<Read, Read>> incoming_pairs;
  std::vector<Read> reads1, reads2;
  while(b->take(incoming_pairs)>=0){
    for(i=0; i<incoming_pairs.size(); ++i){
      reads1.push_back(incoming_pairs[i].first);
      reads2.push_back(incoming_pairs[i].second);
    }
    if(reads1.size()>=VECTOR_SIZE){
#pragma omp parallel for
      for(i=0; i<reads1.size(); ++i){
	std::vector<Read> reads_to_process;
	std::string fastq1, fastq2, fastq_te, read_info;
	reads_to_process.clear();
	reads_to_process.push_back(reads1[i]);
	reads_to_process.push_back(reads2[i]);
	if(process_reads_paired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads_to_process, fastq1, fastq2, fastq_te, read_info)>=0){
	  fastq1_txt_array[omp_get_thread_num()]=fastq1_txt_array[omp_get_thread_num()]+fastq1;
	  fastq2_txt_array[omp_get_thread_num()]=fastq2_txt_array[omp_get_thread_num()]+fastq2;
	  fastq_te_txt_array[omp_get_thread_num()]=fastq_te_txt_array[omp_get_thread_num()]+fastq_te;
	  read_info_txt_array[omp_get_thread_num()]=read_info_txt_array[omp_get_thread_num()]+read_info;
	}
      }
      reads1.clear();
      reads2.clear();
    }
  }
#pragma omp parallel for
  for(i=0; i<reads1.size(); ++i){
    std::vector<Read> reads_to_process;
    std::string fastq1, fastq2, fastq_te, read_info;
    reads_to_process.clear();
    reads_to_process.push_back(reads1[i]);
    reads_to_process.push_back(reads2[i]);
    if(process_reads_paired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads_to_process, fastq1, fastq2, fastq_te, read_info)>=0){
      fastq1_txt_array[omp_get_thread_num()]=fastq1_txt_array[omp_get_thread_num()]+fastq1;
      fastq2_txt_array[omp_get_thread_num()]=fastq2_txt_array[omp_get_thread_num()]+fastq2;
      fastq_te_txt_array[omp_get_thread_num()]=fastq_te_txt_array[omp_get_thread_num()]+fastq_te;
      read_info_txt_array[omp_get_thread_num()]=read_info_txt_array[omp_get_thread_num()]+read_info;
    }
  }


  for(i=0; i<n_threads; ++i){
    rt->fastq1_txt=rt->fastq1_txt+fastq1_txt_array[i];
    rt->fastq2_txt=rt->fastq2_txt+fastq2_txt_array[i];
    rt->fastq_te_txt=rt->fastq_te_txt+fastq_te_txt_array[i];
    rt->read_info_txt=rt->read_info_txt+read_info_txt_array[i];
  }
  return 0;
}



int process_buffer_unpaired(const double &alignment_quality, const double &match_quality, int transposon_length, int host_length, const std::vector<LTR> &TE_references, buffer<Read> *b, std::string* fastq1_txt, std::string* fastq_te_txt, std::string* read_info_txt){
  Read read;
  std::vector<Read> reads, incoming_reads;
  int i, j;
  int n_threads;
#pragma omp parallel
  n_threads=omp_get_num_threads();

  *fastq1_txt="";
  *fastq_te_txt="";
  *read_info_txt="";
  std::string fastq1_txt_array[n_threads];
  std::string fastq_te_txt_array[n_threads];
  std::string read_info_txt_array[n_threads];
  for(i=0; i<n_threads; ++i){
    fastq1_txt_array[i]="";
    fastq_te_txt_array[i]="";
    read_info_txt_array[i]="";
  }
   
  while(b->take(incoming_reads)>=0){
    for(i=0; i<incoming_reads.size(); ++i){
      reads.push_back(incoming_reads[i]);
    }
    if(reads.size()>=VECTOR_SIZE){
#pragma omp parallel for
      for(i=0; i<reads.size(); ++i){
	std::vector<Read> read_to_process;
	std::string fastq1, fastq_te, read_info;
	read_to_process.clear();
	read_to_process.push_back(reads[i]);
	if(process_reads_unpaired(alignment_quality, match_quality, transposon_length, host_length, TE_references, read_to_process, fastq1, fastq_te, read_info)>=0){
	  /*
	    #pragma omp critical
	    {
	    *fastq1_txt=*fastq1_txt+fastq1;
	    *fastq_te_txt=*fastq_te_txt+fastq_te;
	    *read_info_txt=*read_info_txt+read_info;
	    }
	  */
	  fastq1_txt_array[omp_get_thread_num()]=fastq1_txt_array[omp_get_thread_num()]+fastq1;
	  fastq_te_txt_array[omp_get_thread_num()]=fastq_te_txt_array[omp_get_thread_num()]+fastq_te;
	  read_info_txt_array[omp_get_thread_num()]=read_info_txt_array[omp_get_thread_num()]+read_info;
	}
      }
      reads.clear();
    }
  }
#pragma omp parallel for
  for(i=0; i<reads.size(); ++i){
    std::vector<Read> read_to_process;
    std::string fastq1, fastq_te, read_info;
    read_to_process.clear();
    read_to_process.push_back(reads[i]);
    if(process_reads_unpaired(alignment_quality, match_quality, transposon_length, host_length, TE_references, read_to_process, fastq1, fastq_te, read_info)>=0){
      /*
	#pragma omp critical
	{
	*fastq1_txt=*fastq1_txt+fastq1;
	*fastq_te_txt=*fastq_te_txt+fastq_te;
	*read_info_txt=*read_info_txt+read_info;
	}
      */
      fastq1_txt_array[omp_get_thread_num()]=fastq1_txt_array[omp_get_thread_num()]+fastq1;
      fastq_te_txt_array[omp_get_thread_num()]=fastq_te_txt_array[omp_get_thread_num()]+fastq_te;
      read_info_txt_array[omp_get_thread_num()]=read_info_txt_array[omp_get_thread_num()]+read_info;
    }
  }
  for(i=0; i<n_threads; ++i){
    *fastq1_txt=*fastq1_txt+fastq1_txt_array[i];
    *fastq_te_txt=*fastq_te_txt+fastq_te_txt_array[i];
    *read_info_txt=*read_info_txt+read_info_txt_array[i];
  }
  return 0;
}

