// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include "smithwaterman.hpp"
#include "ssw_cpp.h"
#include <iostream>
#include "string.h"
#include <boost/concept_check.hpp>
#include "samdata.hpp"
using namespace std;

void getTrimmedValues(const SmithWatermanResult& smithWatermanResult, const Read& read, int direction, std::string& trimmedRead, std::string& trimmedPhred, int& trimmedPosition, std::string& matchRead, std::string& matchPhred){
  if(direction<0){
    trimmedRead=read.sequence.substr(0, smithWatermanResult.referenceStart);
    trimmedPhred=read.phred.substr(0, smithWatermanResult.referenceStart);
    trimmedPosition=0;
    matchRead=read.sequence.substr(smithWatermanResult.referenceStart);
    matchPhred=read.phred.substr(smithWatermanResult.referenceStart);
  }else{
    trimmedRead=read.sequence.substr(smithWatermanResult.referenceEnd+1);
    trimmedPhred=read.phred.substr(smithWatermanResult.referenceEnd+1);
    trimmedPosition=smithWatermanResult.referenceEnd+1;
    matchRead=read.sequence.substr(0, smithWatermanResult.referenceEnd+1);
    matchPhred=read.phred.substr(0, smithWatermanResult.referenceEnd+1);
  }
}


int smithwaterman(const std::string& shortSequence, const std::string& longSequence, SmithWatermanResult& result){
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;
  aligner.Align(shortSequence.c_str(), longSequence.c_str(), longSequence.size(), filter, &alignment);
  result.score=alignment.sw_score;
  result.relativeScore=((double)result.score)/2/shortSequence.length();
  result.referenceStart=alignment.ref_begin;
  result.referenceEnd=alignment.ref_end;
  result.queryStart=alignment.query_begin;
  result.queryEnd=alignment.query_end;
  result.n_mismatches=alignment.mismatches;
  result.cigarString=alignment.cigar_string;
  return 0; 
}
