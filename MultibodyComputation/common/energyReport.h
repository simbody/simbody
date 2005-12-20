
#ifndef __energyReport_hh__
#define __energyReport_hh__

#include <cdsString.h>

//
// EnergyReport is a simple object that is used to 
// report the results of an energy calculation by a 
// PotentialTerm.  
//
// It holds onto three items--the energy, the deviation,
// and a description string (ie., the name of the potentialTerm
// that generated it)
//

class EnergyReport {
  
public:
  String     name;
  float_type energy;
  float_type deviation;

  EnergyReport(const String     name="",
	       const float_type energy=0.,
	       const float_type deviation=0.) :
    name(name), energy(energy), deviation(deviation) {}

};

#endif /* __energyReport_hh__ */


