#ifdef NOTDEF

#include "avePot.h"
#include <ensembleSimulation.h>
#include <derivList.h>

AvePot::AvePot(rc_Pot      pot,
	       Simulation* simulation)
  : Pot( String("ave_") + pot->potName(), pot->instanceName() ),
    pot(pot)
{
 if ( simulation->type() == String("EnsembleSimulation") ) {
   sim = (EnsembleSimulation*)simulation;
   sharedReports.resize( sim->size() );
 } else
   sim = 0;
} /* constructor */


AvePot::~AvePot()
{} /* destructor */

EnergyReport
AvePot::energyMaybeDerivs(DerivList& derivs,
			  bool       calcDerivs)
{
 int memberIndex = sim->member()->memberIndex();

 //
 // potential terms which are the same for each member
 //
 if ( calcDerivs ) {
   DerivList localDerivs;
   sharedReports[memberIndex] = pot->calcEnergyAndDerivs(localDerivs).energy;
   for (int i=0 ; i<sim->numAtoms() ; i++)
     derivs[sim][i] += scale() * 
		       sim->weight(memberIndex) * 
		       localDerivs[sim->subSim()][i];
 } else
   sharedReports[memberIndex] = pot->calcEnergy().energy;

 // need energy from all members
 sim->barrier(__FILE__,__LINE__); 

 EnergyReport ret(instanceName());
 for (int i=0 ; i<sim->size() ; i++)
   ret.energy += scale() * sim->weight(i) * sharedReports[i];
 
 // don't let shared data get updated before we're finished
 // this should not be necessary
 sim->barrier(__FILE__,__LINE__); 

 return ret;
} /* energyMaybeDerivs */

EnergyReport
AvePot::calcEnergy()
{
 if ( sim ) {
   DerivList tmp;
   return energyMaybeDerivs(tmp,0);
 } else
   return pot->calcEnergy();
} /* calcEnergy */

EnergyReport
AvePot::calcEnergyAndDerivs(DerivList& derivList)
{
 if ( sim )
   return energyMaybeDerivs(derivList,1);
 else 
   return pot->calcEnergyAndDerivs(derivList);
} /* calcEnergyAndDerivs */

#endif
