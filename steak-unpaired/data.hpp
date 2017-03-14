#ifndef DATA_HPP
#define DATA_HPP

#include "samdata.hpp"
#include <string>
#include <vector>


class Read{
public:
  std::string chromosome;
  std::string name;
  int position;
  int quality;
  std::string sequence;
  std::string phred;
  int set(const SAMData& samData);
};

class Pair{
public:
  std::string name;
  Read reads[2];
  int set(const SAMData& samData1, const SAMData& samData2);
  inline Read& operator[](int i){
    return reads[i];
  }
};


#endif
