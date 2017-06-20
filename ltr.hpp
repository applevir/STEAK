#ifndef LTR_HPP
#define LTR_HPP

#include <string>

class LTR{
public:
  std::string name, sequence;
  int direction;
  LTR(std::string i_name, std::string i_sequence, int i_direction){
    name=i_name;
    sequence=i_sequence;
    direction=i_direction;
  }
};

#endif
