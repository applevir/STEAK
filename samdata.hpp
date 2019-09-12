// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


#ifndef SAMDATA_HPP
#define SAMDATA_HPP

#include <string>
int getCigarValue(const std::string& cigarString);

class SAMData{
public:
  std::string name;
  std::string chromosome;
  int position;
  int quality;
  std::string sequence;
  std::string phred;
  std::string cigarString;
  double cigarValue;
  int set(const std::string &line);
};

#endif
