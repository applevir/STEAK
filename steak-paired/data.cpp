#include <math.h>
#include <iostream>
#include "data.hpp"
using namespace std;

int Read::set(const SAMData& samData){
  chromosome=samData.chromosome;
  position=samData.position;
  quality=samData.quality;
  sequence=samData.sequence;
  phred=samData.phred;
  return 0;
}

int Pair::set(const SAMData& samData1, const SAMData& samData2){
  name=samData1.name;
  reads[0].set(samData1);
  reads[1].set(samData2);
  return 0;
}
