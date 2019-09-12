// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "parallelenvironment.hpp"
#include "inputfile.hpp"
#include "outputfile.hpp"
#include <vector>
#include "ltr.hpp"
#include "smithwaterman.hpp"
#include "process.hpp"
#include <boost/circular_buffer.hpp>
#include "buffer.hpp"

#define N_READS_READ_AT_A_TIME 1000
#define BUFFER_SIZE 10000


int main(int argc, char** argv){

  parallel_environment pe(argc, argv);

  namespace bpo=boost::program_options;

  bpo::options_description parameters("Parameters");
  parameters.add_options()
    ("input",             bpo::value<std::string>()->default_value(""), "")
    ("output",            bpo::value<std::string>()->default_value(""), "")
    ("pipe", "")
    ("paired", "")
    ("TE-reference",      bpo::value<std::string>()->default_value(""), "")
    ("alignment-quality", bpo::value<double>()->default_value(.99), "")
    ("match-quality",     bpo::value<double>()->default_value(0.95), "")
    ("transposon-length", bpo::value<int>()->default_value(15), "")
    ("host-length",       bpo::value<int>()->default_value(25), "")
    ("aligned", "")
    ;
  bpo::variables_map variables;
  std::string input, output, te_reference;
  double alignment_quality;
  double match_quality;
  int transposon_length;
  int host_length;

  std::string error_message="\n\
            steak\n\n\
            options [default values]:\n\n\
            --input               Input NGS file  (mandatory without the --pipe option)\n\n\
            --output              Output file (mandatory with the  --pipe option)\n\n\
            --pipe                Takes the input from a pipe. Works only with BAM/SAM and as a single\n\
                                  process.\n\n\
            --paired              Input NGS file(s) are paired-end reads.\n\n\
            --TE-reference        FASTA with sequence of TE/virus of interest.\n\n\
            --alignment-quality   Maximum proportion of ‘M’s in the CIGAR value of a read. Only\n\
            reads of this proportion or lower will be considered.   [.99]\n\n\
            --match-quality       Proportion of the match needed between the TE reference and the read [.95]\n\n\
            --transposon-length   Length of the TE reference that will be searched for within each read. [15]\n\n\
            --host-length         The minimum length that a trimmed read (host flank) can be. Trimmed reads\n\
            smaller than this will be ignored. [25]\n\n\
            --aligned             Uses an aligned genome. Requires a BAM or SAM file. Without this option,\n\
            STEAK will consider that the input consists of a fasta or fastq file and\n\
            the alignment-quality parameter will be ignored.\
            \n\n";


  try{
    bpo::store(bpo::parse_command_line(argc, argv, parameters), variables);
    input=variables["input"].as<std::string>();
    output=variables["output"].as<std::string>();
    te_reference=variables["TE-reference"].as<std::string>();
    alignment_quality=variables["alignment-quality"].as<double>();
    match_quality=variables["match-quality"].as<double>();
    transposon_length=variables["transposon-length"].as<int>();
    host_length=variables["host-length"].as<int>();
  }catch(std::exception& e){
    if(parallel_environment::get_process_num()==0){
      std::cout << error_message;
    }
    return -1;
  }

  bool pipe;
  if(variables.count("pipe")>0){
    pipe=true;
  }else{
    pipe=false;
  }
  bool paired;
  if(variables.count("paired")>0){
    paired=true;
  }else{
    paired=false;
  }
  bool aligned;
  if(variables.count("aligned")>0){
    aligned=true;
  }else{
    aligned=false;
  }

  if(pipe && !aligned){
    std::cout << "Reading data through a pipe requires an aligned file.\n";
  }

  if((pipe && output=="") || (!pipe && input=="") || (pipe && input!="")){
    if(parallel_environment::get_process_num()==0){
      std::cout << error_message;
    }
    return -1;
  }

  if(te_reference==""){
    if(parallel_environment::get_process_num()==0){
      std::cout << "A TE reference file is required.\n";
    }
    return -1;
  }


  /*
    std::ifstream inBait(te_reference);
    if(!inBait.is_open()){
    std::cout << "Unable to open the TE reference file\n";
    return -1;
    }
  */

  fasta_file* TE_references_file;
  try{
    TE_references_file=new fasta_file(te_reference, true);
  }catch(std::exception e){
    if(parallel_environment::get_process_num()==0){
      std::cout << "Unable to open the TE reference file.\n";
    }
    return -1;
  }
  std::string line, bait_left, bait_right, bait_left_complement, bait_right_complement;
  std::vector<LTR> TE_references;
  Read reference_read;


  while(TE_references_file->get_read(reference_read)){
    if(reference_read.sequence.length()<transposon_length){
      if(parallel_environment::get_process_num()==0){
	std::cout << "The TE reference is not long enough.\n";
      }
      return -1;
    }
    bait_left=reference_read.sequence.substr(0, transposon_length);
    bait_right=reference_read.sequence.substr(reference_read.sequence.length()-transposon_length);
    bait_left_complement=create_complement(bait_left);
    bait_right_complement=create_complement(bait_right);
    TE_references.push_back(LTR(reference_read.name, bait_left, -1));
    TE_references.push_back(LTR(reference_read.name, bait_right, 1));
    TE_references.push_back(LTR(reference_read.name, bait_left_complement, 1));
    TE_references.push_back(LTR(reference_read.name, bait_right_complement, -1));
  }



  std::string name;
  name=input.substr(0, input.rfind("."));
  if(output==""){
    output=name;
  }else{
    name=output;
  }
  std::string fastq1, fastq2, fastq_te, read_info;
  std::string fastq1_txt, fastq2_txt, fastq_te_txt, read_info_txt;


  if(aligned){
 
    if(pipe){

      if(parallel_environment::get_num_processes()>1){
	std::cout << "The reading of a BAM file through a pipe is allowed with only one process\n";
      }

      if(paired){



	buffer<std::pair<Read, Read>> b(BUFFER_SIZE);
	bam_pipe sf;
	Read read1, read2;
	std::vector<std::pair<Read, Read>> pairs_accumulated;
	result_texts rt;

	boost::thread t(&process_buffer_paired, alignment_quality, match_quality, transposon_length, host_length, TE_references, &b, &rt);



	while(sf.get_read(read1) && sf.get_read(read2)){
	    pairs_accumulated.push_back(std::pair<Read, Read>(read1, read2));
            if(pairs_accumulated.size()>=N_READS_READ_AT_A_TIME/2){
              b.put(pairs_accumulated);
              pairs_accumulated.clear();
	    }
	}
        b.put(pairs_accumulated);
	b.stop_waiting();
	t.join();

	output_file out_fastq1(name+".1.fastq", rt.fastq1_txt.length());
	out_fastq1.write(rt.fastq1_txt);
	output_file out_fastq2(name+".2.fastq", rt.fastq2_txt.length());
	out_fastq2.write(rt.fastq2_txt);
	output_file out_fastq_te(name+".te.fastq", rt.fastq_te_txt.length());
	out_fastq_te.write(rt.fastq_te_txt);
	output_file out_read_info(name+".dat", rt.read_info_txt.length());
	out_read_info.write(rt.read_info_txt);
	return 0;
      }else{

	buffer<Read> b(BUFFER_SIZE);
	bam_pipe sf;
	Read read;
	std::vector<Read> reads_accumulated;
	fastq1_txt="";
	fastq_te_txt=""; 
	read_info_txt="";
	boost::thread t{&process_buffer_unpaired, alignment_quality, match_quality, transposon_length, host_length, TE_references, &b, &fastq1_txt, &fastq_te_txt, &read_info_txt};


	while(sf.get_read(read)){
          reads_accumulated.push_back(read);
          if(reads_accumulated.size()>=N_READS_READ_AT_A_TIME){
	    b.put(reads_accumulated);
            reads_accumulated.clear();
	  }
	}
	b.put(reads_accumulated);
	b.stop_waiting();
	t.join();
	output_file out_fastq1(name+".fastq", fastq1_txt.length());
	out_fastq1.write(fastq1_txt);
	output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
	out_fastq_te.write(fastq_te_txt);
	output_file out_read_info(name+".dat", read_info_txt.length());
	out_read_info.write(read_info_txt);
	return 0;
      }
      /*

	if(paired){
	bam_pipe sf;
	Read read1, read2;
	std::vector<Read> reads;
	bool first_time=true;
	while(sf.get_read(read1)){
	sf.get_read(read2);
	if(first_time && read1.name.substr(0, read1.name.rfind("/"))!=read2.name.substr(0, read2.name.rfind("/"))){
	read1=read2;
	sf.get_read(read2);
	}
	reads.clear();
	if(read1.name.substr(0, read1.name.rfind("/"))!=read2.name.substr(0, read2.name.rfind("/"))){
	std::cout << " Paired-ends files must be sorted by name.\n";
	return -1;
	}
	reads.push_back(read1);
	reads.push_back(read2);
	if(process_reads_paired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads, fastq1, fastq2, fastq_te, read_info)>=0){
	fastq1_txt=fastq1_txt+fastq1;
	fastq2_txt=fastq2_txt+fastq2;
	fastq_te_txt= fastq_te_txt+fastq_te;
	read_info_txt=read_info_txt+read_info;
	}
	}
	output_file out_fastq1(name+".1.fastq", fastq1_txt.length());
	out_fastq1.write(fastq1_txt);
	output_file out_fastq2(name+".2.fastq", fastq2_txt.length());
	out_fastq2.write(fastq2_txt);
	output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
	out_fastq_te.write(fastq_te_txt);
	output_file out_read_info(name+".dat", read_info_txt.length());
	out_read_info.write(read_info_txt);
	}else{
	bam_pipe sf;
	Read read1;
	std::vector<Read> reads;
	bool first_time=true;
	while(sf.get_read(read1)){
	reads.clear();
	reads.push_back(read1);
	if(process_reads_unpaired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads, fastq1, fastq_te, read_info)>=0){
	fastq1_txt=fastq1_txt+fastq1;
	fastq_te_txt= fastq_te_txt+fastq_te;
	read_info_txt=read_info_txt+read_info;
	}
	}
	output_file out_fastq1(name+".fastq", fastq1_txt.length());
	out_fastq1.write(fastq1_txt);
	output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
	out_fastq_te.write(fastq_te_txt);
	output_file out_read_info(name+".dat", read_info_txt.length());
	out_read_info.write(read_info_txt);

	}
      */


    }else{
      if(paired){
	sam_file sf(input);
	Read read1, read2;
	std::vector<Read> reads;
	bool first_time=true;
	while(sf.get_read(read1)){
	  sf.get_read(read2);
	  if(first_time && read1.name.substr(0, read1.name.rfind("/"))!=read2.name.substr(0, read2.name.rfind("/"))){
	    read1=read2;
	    sf.get_read(read2);
	  }
	  reads.clear();
	  if(read1.name.substr(0, read1.name.rfind("/"))!=read2.name.substr(0, read2.name.rfind("/"))){
	    std::cout << " Paired-ends files must be sorted by name.\n";
	    return -1;
	  }
	  reads.push_back(read1);
	  reads.push_back(read2);
	  if(process_reads_paired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads, fastq1, fastq2, fastq_te, read_info)>=0){
	    fastq1_txt=fastq1_txt+fastq1;
	    fastq2_txt=fastq2_txt+fastq2;
	    fastq_te_txt= fastq_te_txt+fastq_te;
	    read_info_txt=read_info_txt+read_info;
	  }
	}
	output_file out_fastq1(name+".1.fastq", fastq1_txt.length());
	out_fastq1.write(fastq1_txt);
	output_file out_fastq2(name+".2.fastq", fastq2_txt.length());
	out_fastq2.write(fastq2_txt);
	output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
	out_fastq_te.write(fastq_te_txt);
	output_file out_read_info(name+".dat", read_info_txt.length());
	out_read_info.write(read_info_txt);
	return 0;
      }else{
	sam_file sf(input);
	Read read1;
	std::vector<Read> reads;
	bool first_time=true;
	while(sf.get_read(read1)){
	  reads.clear();
	  reads.push_back(read1);
	  if(process_reads_unpaired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads, fastq1, fastq_te, read_info)>=0){
	    fastq1_txt=fastq1_txt+fastq1;
	    fastq_te_txt= fastq_te_txt+fastq_te;
	    read_info_txt=read_info_txt+read_info;
	  }
	}
	output_file out_fastq1(name+".fastq", fastq1_txt.length());
	out_fastq1.write(fastq1_txt);
	output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
	out_fastq_te.write(fastq_te_txt);
	output_file out_read_info(name+".dat", read_info_txt.length());
	out_read_info.write(read_info_txt);
	return 0;
      }


    }
  }else{

    if(paired){
      fastq_file sf(input);
      Read read1, read2;
      std::vector<Read> reads;
      sf.get_start_position();
      sf.get_pair_beginning();
      while(sf.get_read(read1) && sf.get_read(read2, true)){
	
	if(read1.name.substr(0, read1.name.rfind("/"))!=read2.name.substr(0, read2.name.rfind("/"))){
	  std::cout << "Paired-ends files must be sorted by name.\n";
	  return -1;
	}
	reads.clear();
	reads.push_back(read1);
	reads.push_back(read2);
	if(process_reads_paired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads, fastq1, fastq2, fastq_te, read_info)>=0){
	  fastq1_txt=fastq1_txt+fastq1;
	  fastq2_txt=fastq2_txt+fastq2;
	  fastq_te_txt= fastq_te_txt+fastq_te;
	  read_info_txt=read_info_txt+read_info;
	}
      }
      output_file out_fastq1(name+".1.fastq", fastq1_txt.length());
      out_fastq1.write(fastq1_txt);
      output_file out_fastq2(name+".2.fastq", fastq2_txt.length());
      out_fastq2.write(fastq2_txt);
      output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
      out_fastq_te.write(fastq_te_txt);
      output_file out_read_info(name+".dat", read_info_txt.length());
      out_read_info.write(read_info_txt);
      return 0;
    }else{
      fastq_file sf(input);
      Read read1;
      std::vector<Read> reads;
      sf.get_start_position();
      sf.get_pair_beginning();
      while(sf.get_read(read1)){
	reads.clear();
	reads.push_back(read1);
	if(process_reads_unpaired(alignment_quality, match_quality, transposon_length, host_length, TE_references, reads, fastq1, fastq_te, read_info)>=0){
	  fastq1_txt=fastq1_txt+fastq1;
	  fastq_te_txt= fastq_te_txt+fastq_te;
	  read_info_txt=read_info_txt+read_info;
	}
      }
      output_file out_fastq1(name+".fastq", fastq1_txt.length());
      out_fastq1.write(fastq1_txt);
      output_file out_fastq_te(name+".te.fastq", fastq_te_txt.length());
      out_fastq_te.write(fastq_te_txt);
      output_file out_read_info(name+".dat", read_info_txt.length());
      out_read_info.write(read_info_txt);
      return 0;
    }


  }




  return 0;

}
