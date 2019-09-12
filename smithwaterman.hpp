// STEAK   Specific Transposable Element Aligner (HERV-K) 
// Copyright (C) 2015 Cindy Santander, Philippe Gambron
// Distributed under the GNU General Public License version 3 (GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)


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

void getTrimmedValues(const SmithWatermanResult& smithWatermanResult, const Read& read, int direction, std::string& trimmedRead, std::string& trimmedPhred, int& trimmedPosition, std::string& matchRead, std::string& matchPhred);

int smithwaterman(const std::string& shortSequence, const std::string& longSequence, SmithWatermanResult& result);

#endif
