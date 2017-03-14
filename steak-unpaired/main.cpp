#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <string>
#include "samdata.hpp"
#include <vector>
#include <omp.h>
#include "smithwaterman.hpp"
#include "data.hpp"
#include "ltr.hpp"
//#include "checkwithsubstrings.hpp"
//#include "checkdata.hpp"
#include <map>
using namespace std;

string create_complement(const string& sequence){
  string s, c;
  map<string, string> corr;
  corr["A"]="T";
  corr["G"]="C";
  corr["C"]="G";
  corr["T"]="A";
  for(int i=sequence.length()-1; i>=0; --i){
    s=s+corr[sequence.substr(i, 1)];
  }
  return s;
}


int main(int argc, char** argv){

  int i, j, i_max, i_read, i_LTR, bait_length;
  if(argc<7){
    return -1;
  }
  string name, fileName(*(argv+1)), baitFileName(*(argv+5));  
  ifstream in(fileName);
  istream *input;
  if(in.is_open()){
    input=&in;
    name=fileName.substr(0, fileName.rfind("."));
  }else{
    input=&cin;
    name=fileName;
  }
  
  string fastqName, dataName, dataNameOrig;
  fastqName=name+".fastq";
  dataName=name+".txt";
  dataNameOrig=name+".orig.txt";
  ofstream fastqOut(fastqName);
  ofstream dataOut(dataName);
  ofstream dataOutOrig(dataNameOrig);
  
  double cigarCeiling=boost::lexical_cast<double>(*(argv+2));
  double smithWatermanThreshold=boost::lexical_cast<double>(*(argv+3));
  int minimalTrimmedLength=boost::lexical_cast<int>(*(argv+4));
//  int subStringLength=boost::lexical_cast<int>(*(argv+5));
//  double subStringThreshold=boost::lexical_cast<double>(*(argv+6));
//  int minimalLengthIfAtTheEnd=boost::lexical_cast<int>(*(argv+5));
  bait_length=boost::lexical_cast<double>(*(argv+6));
  
  dataOut << "# Cigar value ceiling:      " << cigarCeiling << "\n";
  dataOut << "# Smith-Waterman threshold: " << smithWatermanThreshold << "\n";
  dataOut << "# Minimal trimmed length:   " << minimalTrimmedLength << "\n";
  //  dataOut << "# Minimal length if at extremity: " << minimalLengthIfAtTheEnd << "\n\n";

  string line, bait_left, bait_right, bait_left_complement, bait_right_complement;
  ifstream inBait(baitFileName);
  if(!inBait.is_open()){
    return -1;
  }
  vector<LTR> LTRs;
  while(getline(inBait, line)){
    if(line.length()>0 && line.substr(0, 1)!=">"){
      bait_left=line.substr(0, bait_length);
      bait_right=line.substr(line.length()-bait_length);
      bait_left_complement=create_complement(bait_left);
      bait_right_complement=create_complement(bait_right);
      LTRs.push_back(LTR(bait_left, -1));
      LTRs.push_back(LTR(bait_right, 1));
      LTRs.push_back(LTR(bait_left_complement, 1));
      LTRs.push_back(LTR(bait_right_complement, -1));
    }
  }
  inBait.close();
  
  
  //  vector<LTR> LTRs;
  // LTRs.push_back(LTR("TGTGGGGAAAAGCAAGAGAG", -1));
  //LTRs.push_back(LTR("AGGGGCAACCCACCCCTACA", 1));
  //LTRs.push_back(LTR("TGTAGGGGTGGGTTGCCCCT", -1));
  //LTRs.push_back(LTR("CTCTCTTGCTTTTCCCCACA", 1));
  
  SmithWatermanResult smithWatermanResult[LTRs.size()];
  SAMData samData;
  double maxScore;
  // double subStringProportion;
  
  string  trimmedRead[LTRs.size()], trimmedQuality[LTRs.size()];
  int trimmedPosition[LTRs.size()];
  double score[LTRs.size()];
  Read read;
  int i_maxScore;  
  while(getline(*input, line)){
    samData.set(line);
    read.set(samData);
    if(samData.cigarValue<=cigarCeiling){
#pragma omp parallel for
      for(i=0; i<LTRs.size(); ++i){
	smithwaterman(LTRs[i].sequence, samData.sequence, smithWatermanResult[i]);
	getTrimmedValues(smithWatermanResult[i], read, LTRs[i].direction, trimmedRead[i], trimmedQuality[i], trimmedPosition[i]);
	score[i]=smithWatermanResult[i].relativeScore;
      }
      maxScore=0;
      i_maxScore=-1;
      for(i=0; i<LTRs.size(); ++i){
        if(score[i]>=smithWatermanThreshold && trimmedRead[i].length()>=minimalTrimmedLength){
          if(score[i]>maxScore){
            i_maxScore=i;
            maxScore=score[i];
          }
        }
      }
            
      if(i_maxScore>=0){

	fastqOut << "@" << read.name << "\n";
	fastqOut << trimmedRead[i_maxScore] << "\n";
	fastqOut << "+\n";
	fastqOut << trimmedQuality[i_maxScore] << "\n";

	dataOutOrig << read.name << " " << read.chromosome << " " << read.position << " " << read.sequence << "\n";

	dataOut << read.name  << " " << read.chromosome << " " << read.position+trimmedPosition[i_maxScore] << " " << trimmedRead[i_maxScore] << " " << maxScore << "\n";
      }
    }
   
  }
  
  fastqOut.close();
  dataOut.close();
  dataOutOrig.close();
  return 0; 
}
