// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include "inputfile.hpp"
#include "parallelenvironment.hpp"
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include "data.hpp"

class file_exception: public std::exception{
  virtual const char* what() const throw(){
    return "Unable to open file";
  }
} file_exception;


fasta_file::fasta_file(const std::string filename, bool sequential){
  unsigned long long size;
  in.open(filename, std::ios::in);
  if(!in.is_open()){
    throw file_exception;
  }
  int pos=filename.rfind(".");
  if(parallel_environment::get_process_num()==0 && (pos==std::string::npos || pos==filename.length()-1 || !(filename.substr(pos+1)=="fasta" || filename.substr(pos+1)=="fa" || filename.substr(pos+1)=="fastq" || filename.substr(pos+1)=="fq"))){
    std::cout << "The TE reference file must be a fasta file.\n";
  }
  in.seekg(0, in.end);
  size=in.tellg();
  int num_processes, process_num;

  if(sequential){
    num_processes=1;
    process_num=0;
  }else{
    num_processes=parallel_environment::get_num_processes();
    process_num=parallel_environment::get_process_num();
  }
  local_size=ceil(((double)size)/num_processes);
  local_beginning=local_size*process_num;
  local_end=local_beginning+local_size;
  in.seekg(local_beginning);
  local_position=local_beginning;
}

fasta_file::~fasta_file(){
  in.close();
}

bool fasta_file::get_read(Read& read){

  bool keep=false, within_sequence;
  std::string line;
  std::string sequence="", name;
  while((local_position<local_end || keep) && std::getline(in, line)){
    boost::trim(line);
    if(keep && line.length()>0 && (line.substr(0,1)==">" || line.substr(0,1)=="@")){
      in.seekg(local_position);
      read.set(name, sequence);
      return true;
    }
    if(line.length()>0 && line.substr(0,1)=="+"){
      within_sequence=false;
    }
    if(keep && within_sequence){
      sequence=sequence+line;
    }
    if(!keep && line.length()>0 && (line.substr(0,1)==">" || line.substr(0,1)=="@")){
      keep=true;
      within_sequence=true;
      name=line.substr(1);
      boost::trim(name);
    }
    local_position=in.tellg();
  }
  if(in.eof() && keep){
    read.set(name, sequence);
    return true;
  }
  return false;
}



fastq_file::fastq_file(const std::string& filename){
  unsigned long long size;
  in.open(filename, std::ios::in);
  if(!in.is_open()){
    throw file_exception;
  }
  int pos=filename.rfind(".");
  if(parallel_environment::get_process_num()==0 && (pos==std::string::npos || pos==filename.length()-1 || !(filename.substr(pos+1)=="fastq" || filename.substr(pos+1)=="fq"))){
    std::cout << "The input file must be a fastq file.\n";
  }
  in.seekg(0, in.end);
  size=in.tellg();
  local_size=ceil(((double)size)/parallel_environment::get_num_processes());
  local_beginning=local_size*parallel_environment::get_process_num();
  local_end=local_beginning+local_size;
  in.seekg(local_beginning);
  local_position=local_beginning;
  std::string line;
}

fastq_file::~fastq_file(){
  in.close();
}

bool fastq_file::get_start_position(bool recursive){
  if(parallel_environment::get_process_num()>0){
    std::string line;
    std::vector<unsigned long long> positions;
    std::vector<std::string> lines;
    do{
      positions.push_back(in.tellg());
      if(!std::getline(in, line)){
        return false;
      }
    }while(line!="+");
    if(positions.size()>3){
      in.seekg(positions[positions.size()-3]);   
    }else{
        if(recursive){
          get_start_position(false);
        }
    }
  }
  return true;
}

bool fastq_file::get_pair_beginning(){
    unsigned long long pos1, pos2;
    Read read1, read2;
    pos1=in.tellg();
    if(!get_read(read1)){
      return false;
    }
    pos2=in.tellg();
    if(!get_read(read2)){
      return false;
    }
    if(read1.name.substr(0, read1.name.rfind("/"))==read2.name.substr(0, read2.name.rfind("/"))){
	  in.seekg(pos1);
	}else{
      in.seekg(pos2);  
    }
    return true;
}


bool fastq_file::get_read(Read& read, bool second_read){
  if(local_position>=local_end+1 &!second_read){
	  return false;
  }
  int i_line=0;
  std::string name, sequence, quality, line;
  while(std::getline(in, line)){
	switch(i_line){
	    case 0:
	      	name=line.substr(1);
	      	break;
	    case 1:
	    	sequence=line;
           	break;
	    case 3:
	    	quality=line;
	    	read.set(name, sequence, quality);
	    	local_position=in.tellg();
           	return true;
	    default:
	    	break;
	    }
		++i_line;
	}
	return false;
  }



sam_file::sam_file(const std::string& filename){
	  unsigned long long size;
	  in.open(filename, std::ios::in);
	  if(!in.is_open()){
	    throw file_exception;
	  }
	  int pos=filename.rfind(".");
	  if(parallel_environment::get_process_num()==0 && (pos==std::string::npos || pos==filename.length()-1 || !(filename.substr(pos+1)=="sam"))){
	    std::cout << "The input file must be a SAM file. BAM files can be fed through samtools view and a pipe\n";
	  }
	  in.seekg(0, in.end);
	  size=in.tellg();
	  local_size=ceil(((double)size)/parallel_environment::get_num_processes());
	  local_beginning=local_size*parallel_environment::get_process_num();
	  local_end=local_beginning+local_size;
	  in.seekg(local_beginning);
	  local_position=local_beginning;
          std::string line;
          if(parallel_environment::get_process_num()>0){
		getline(in, line);
          }else{
              unsigned long long current_position;
              do{
                  current_position=in.tellg();
                  getline(in, line);
              }while(line.substr(0, 1)=="@");
              in.seekg(current_position);
          }

}

sam_file::~sam_file(){
     in.close();
}

bool sam_file::get_read(Read& read){
  std::string line;
  if(!std::getline(in, line) || local_position>=local_end){
	  return false;
  }
  SAMData samData;
  samData.set(line);
  read.set(samData);
  local_position=in.tellg(); //local_position+line.length()+1;
  return true;
}

bam_pipe::bam_pipe(){
  input=&std::cin;	
}

int bam_pipe::get_read(Read& read){
  std::string line;
  if(!std::getline(*input, line)){
    return false;
  }
  SAMData samData;
  samData.set(line);
  read.set(samData);
  return true;
}
