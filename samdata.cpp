// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#include "samdata.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
using namespace std;

int getCigarValue(const string& cigarString){
  if(cigarString=="*"){
    return 0;
  }
  boost::char_separator<char> separators("", "MIDNSHP=X");
  boost::tokenizer<boost::char_separator<char>> tokens(cigarString, separators);
  boost::tokenizer<boost::char_separator<char>>::iterator it;
  int value, total=0;
  string code;
  for(it=tokens.begin(); it!=tokens.end(); ++it){
    value=boost::lexical_cast<int>(*it);
    ++it;
    code=*it;
    if(code=="M"){
      total=total+value;
    }
  }
  return total;
}

int SAMData::set(const string& line){
  vector<string> elements;
  boost::split(elements, line, boost::is_any_of("\t"));
  name=elements[0]; 
  chromosome=elements[2];
  position=boost::lexical_cast<int>(elements[3]);
  quality=boost::lexical_cast<int>(elements[4]);
  sequence=elements[9];
  phred=elements[10];
  cigarString=elements[5];
  cigarValue=((double)getCigarValue(cigarString))/sequence.length();
  return 0;
}
