
#ifndef __xplorPot_hh__
#define __xplorPot_hh__

#include <pot.h>

//
// class which wraps an xplor potential term
//

class Simulation;
class XplorSimulation;

class XplorPot : public Pot {
public:
  XplorPot(const CDSString&     instanceName,
		   Simulation* defaultSimulation);
  
  virtual EnergyReport calcEnergy();
  virtual EnergyReport calcEnergyAndDerivs(DerivList&);
private:
  int potIndex; // index into xplor qener and nener arrays
  XplorSimulation* sim;

  XplorPot();                        //inaccessible
  XplorPot(const XplorPot&);
  void operator=(const XplorPot&);
};

#endif /* __xplorPot_hh__ */

