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
