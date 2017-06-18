#include <math.h>
#include <iostream>
#include "data.hpp"
using namespace std;

int Read::set(const SAMData& samData){
  name=samData.name;
  chromosome=samData.chromosome;
  position=samData.position;
  quality=samData.quality;
  cigarValue=samData.cigarValue;
  cigarString=samData.cigarString;
  sequence=samData.sequence;
  phred=samData.phred;
  return 0;
}
/*
int Read::set(const std::string& i_chromosome, const std::string& i_name, const int i_position, const int i_quality, const std::string& i_sequence, const std::string& i_phred){
	chromosome=i_chromosome;
	name=i_name;
	position=i_position;
	quality=i_quality;
	sequence=i_sequence;
	phred=i_phred;
	return 0;
}
*/
int Read::set(const std::string& i_name, const std::string& i_sequence, const std::string& i_phred){
	chromosome="";
	name=i_name;
	position=0;
	quality=0;
        cigarString="";
	sequence=i_sequence;
	phred=i_phred;
	return 0;
}
/*
int Pair::set(const SAMData& samData1, const SAMData& samData2){
  name=samData1.name;
  reads[0].set(samData1);
  reads[1].set(samData2);
  return 0;
}
*/
