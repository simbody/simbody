
#include "derivList.h"

#include <simulation.h>

DerivList::DerivList()
{
}

void
DerivList::init(const Simulation* sim)
{
 if ( derivList.size() <= sim->id() )
   derivList.resize( sim->id()+1 );

 if ( derivList[sim->id()].size() != sim->numAtoms() )
   derivList[sim->id()].resize( sim->numAtoms() );

 static Vec3 zeroVec( (float_type)0. );
 for (int i=0 ; i<sim->numAtoms() ; i++)
   derivList[sim->id()](i) = zeroVec;
} /* init */


void
DerivList::clear()
{
 Vec3 zero(0.);
 for (int i=0 ; i<derivList.size() ; i++)
   for (int j=0 ; j<derivList[i].size() ; j++)
     derivList[i][j] = zero;
} /* clear */
