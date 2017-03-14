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
  
  string dataName;
  string fastqName1, fastqName2;
  fastqName1=name+".1.fastq";
  fastqName2=name+".2.fastq";
  dataName=name+".txt";
  ofstream fastqOut1(fastqName1);
  ofstream fastqOut2(fastqName2);
  ofstream dataOut(dataName);
  
  double cigarCeiling=boost::lexical_cast<double>(*(argv+2));
  double smithWatermanThreshold=boost::lexical_cast<double>(*(argv+3));
  int minimalTrimmedLength=boost::lexical_cast<int>(*(argv+4));
  //  int subStringLength=boost::lexical_cast<int>(*(argv+5));
  //  double subStringThreshold=boost::lexical_cast<double>(*(argv+6));
  bait_length=boost::lexical_cast<double>(*(argv+6));
  dataOut << "# Cigar value ceiling:      " << cigarCeiling << "\n";
  dataOut << "# Smith-Waterman threshold: " << smithWatermanThreshold << "\n";
  dataOut << "# Minimal trimmed length:   " << minimalTrimmedLength << "\n\n";
  // dataOut << "# Substring length:   " << subStringLength << "\n\n";
  // dataOut << "# Substring threshold:   " << subStringThreshold << "\n\n";
  
  

  ifstream inBait(baitFileName);
  if(!inBait.is_open()){
    return -1;
  }
  string line, bait_left, bait_right, bait_left_complement, bait_right_complement;
 
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
  
  // for(i=0; i<LTRs.size(); ++i)
  // std::cout << LTRs[i].sequence << " " << LTRs[i].direction << "\n";  
  // return 0;
 // LTRs.push_back(LTR("TGTGGGGAAAAGCAAGAGAG", -1));
 // LTRs.push_back(LTR("AGGGGCAACCCACCCCTACA", 1));
 // LTRs.push_back(LTR("TGTAGGGGTGGGTTGCCCCT", -1));
 // LTRs.push_back(LTR("CTCTCTTGCTTTTCCCCACA", 1));

  
  SmithWatermanResult smithWatermanResult[2*LTRs.size()];
  SAMData samData1, samData2;
  double score, maxScore, subStringProportion;
  Pair pair;
  string line1, line2, trimmedRead[2*LTRs.size()], trimmedQuality[2*LTRs.size()];
  int trimmedPosition[2*LTRs.size()];
  double absoluteScore[2*LTRs.size()];
  int i_pair=0;
  int i_max_for_read[2];
  int i_for_file1, i_for_file2, i_other_read;
  int printed_read;
  double max_score_for_read[2];
  while(getline(*input, line1) && getline(*input, line2)){
    samData1.set(line1);
    samData2.set(line2);
    pair.set(samData1, samData2);
    if(samData1.cigarValue<=cigarCeiling || samData2.cigarValue<=cigarCeiling){
#pragma omp parallel for
      for(i=0; i<LTRs.size()*2; ++i){
	  smithwaterman(LTRs[i/2].sequence, pair[i%2].sequence, smithWatermanResult[i]);
	  getTrimmedValues(smithWatermanResult[i], pair[i%2], LTRs[i/2].direction, trimmedRead[i], trimmedQuality[i], trimmedPosition[i]);
	  absoluteScore[i]=smithWatermanResult[i].score/2;
      }

      

      i_max_for_read[0]=-1;
      i_max_for_read[1]=-1;
      max_score_for_read[0]=0;
      max_score_for_read[1]=0;
      
      for(i=0; i<LTRs.size()*2; ++i){
	i_read=i%2;
	i_LTR=i/2;
	
	if(trimmedRead[i].length()>=minimalTrimmedLength){
	  if(smithWatermanResult[i].relativeScore>=smithWatermanThreshold){
	    if(absoluteScore[i]>max_score_for_read[i_read]){
	      i_max_for_read[i_read]=i;
	      max_score_for_read[i_read]=absoluteScore[i];
	    }
	  }
	}
	
      }
      
      if(i_max_for_read[0]>=0 || i_max_for_read[1]>=0){
	
	if(i_max_for_read[0]>=0){
          i_for_file1=i_max_for_read[0];
          i_for_file2=i_max_for_read[1];
          i_other_read=1;
	  printed_read=0;
        }else{
          i_for_file1=i_max_for_read[1];
          i_for_file2=i_max_for_read[0];
          i_other_read=0;
	  printed_read=1;
        }
		
	fastqOut1 << "@" << pair.name << "\n";
	fastqOut1 << trimmedRead[i_for_file1] << "\n";
	fastqOut1 << "+\n";
	fastqOut1 << trimmedQuality[i_for_file1] << "\n";
	
	if(i_for_file2>=0){
	  fastqOut2 << "@" << pair.name << "\n";
	  fastqOut2 << trimmedRead[i_for_file2] << "\n";
	  fastqOut2 << "+\n";
	  fastqOut2 << trimmedQuality[i_for_file2] << "\n";
	}else{
	  fastqOut2 << "@" << pair.name << "\n";
	  fastqOut2 << pair[i_other_read].sequence << "\n";
	  fastqOut2 << "+\n";
	  fastqOut2 << pair[i_other_read].phred << "\n";	   
	}
	
      
	dataOut << pair.name  << " " << i_read << " " << pair[printed_read].chromosome << " " << pair[printed_read].position+trimmedPosition[i_for_file1] << " " << trimmedRead[i_for_file1] << "   " << pair[0].chromosome << " " << pair[0].position << " " << pair[0].sequence << " " << pair[1].chromosome << " " << pair[1].position << " " << pair[1].sequence << "\n";
	
      }
      
    }
    ++i_pair;
    if((i_pair % 10000) ==0){
      cout << i_pair << "\n";
    }
  }
  
  fastqOut1.close();
  fastqOut2.close();
  dataOut.close();
  return 0; 
}
