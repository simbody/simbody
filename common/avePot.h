#ifdef NOTDEF

#ifndef __avePot_hh__
#define __avePot_hh__

#include <pot.h>

#include <ensembleSimulation.h>
class DerivList;

class AvePot : public Pot {
public:
  AvePot(rc_Pot      pot,
	 Simulation* defaultSimulation);

  ~AvePot();

  virtual EnergyReport calcEnergy();
  virtual EnergyReport calcEnergyAndDerivs(DerivList&);

private:
  EnsembleSimulation* sim;
  rc_Pot pot;
  CDSVector<float_type,0,EnsembleSimulation::SharedAlloc> sharedReports;


  EnergyReport energyMaybeDerivs(DerivList& derivs,
				 bool       calcDerivs);

  AvePot();                        //inaccessible
  AvePot(const AvePot&);
  void operator=(const AvePot&);
};

#endif /* __avePot_hh__ */

#endif

