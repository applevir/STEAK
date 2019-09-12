// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


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
