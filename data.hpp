// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#ifndef DATA_HPP
#define DATA_HPP

#include "samdata.hpp"
#include <string>
#include <vector>
#include <map>

class Read{
public:
  std::string chromosome;
  std::string name;
  int position;
  int quality;
  double cigarValue;
  std::string sequence;
  std::string phred;
  std::string cigarString;
  int set(const SAMData& samData);
  //int set(const std::string& i_chromosome, const std::string& i_name, const int i_position, const int i_quality, const std::string& i_sequence, const std::string& i_phred);
  int set(const std::string& i_name, const std::string& i_sequence, const std::string& i_phred="");
};
/*
class Pair{
public:
  std::string name;
  Read reads[2];
  int set(const SAMData& samData1, const SAMData& samData2);
  inline Read& operator[](int i){
    return reads[i];
  }
};
*/
inline std::string create_complement(const std::string& sequence){
  std::string s, c;
  std::map<std::string, std::string> corr;
  corr["A"]="T";
  corr["G"]="C";
  corr["C"]="G";
  corr["T"]="A";
  for(int i=sequence.length()-1; i>=0; --i){
    s=s+corr[sequence.substr(i, 1)];
  }
  return s;
}

#endif
