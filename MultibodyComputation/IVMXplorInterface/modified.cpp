
#include "modified.h"
#include <simulation.h>
#include <cdsExcept.h>

void
ModifiedBase::registerTo(const Simulation* sim)
{
 sim->addDependent(this);
 registeredSimulations.append(sim);
} /* registerTo */

void
ModifiedBase::unRegister(const Simulation* sim)
{
 if (  !registeredSimulations.contains(sim) )
   throw CDS::exception(CDSString("ModifiedBase::unRegister: ") +
			"mismatch in registration status");

 registeredSimulations.remove( registeredSimulations.getIndex(sim) );
} /* unRegister */


ModifiedBase::~ModifiedBase()
{
 for (int i=0 ; i<registeredSimulations.size() ; i++)
   if ( registeredSimulations[i] )
     registeredSimulations[i]->removeDependent(this);
} /* destructor */

