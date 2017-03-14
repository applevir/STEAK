#ifndef SMITHWATERMAN_HPP
#define SMITHWATERMAN_HPP

#include <string>
#include "data.hpp"


class SmithWatermanResult{
public:
  int score;
  double relativeScore;
  int queryStart, queryEnd;
  int referenceStart, referenceEnd;
  int n_mismatches;
  std::string cigarString;
};

void getTrimmedValues(const SmithWatermanResult& smithWatermanResult, const Read& read, int direction, std::string& trimmedRead, std::string& trimmedQuality, int& trimmedPosition);

int smithwaterman(const std::string& shortSequence, const std::string& longSequence, SmithWatermanResult& result);

#endif