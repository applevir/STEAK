#ifndef LTR_HPP
#define LTR_HPP

#include <string>

class LTR{
public:
  std::string sequence;
  int direction;
  LTR(std::string i_sequence, int i_direction){
    sequence=i_sequence;
    direction=i_direction;
  }
};

#endif