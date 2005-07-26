
#ifdef NOTDEF

#include "xplorPot.h"
#include <sthead.h>
#include <cdsExcept.h>
#include <xplorVars.h>
#include <derivList.h>
#include <cdsAuto_arr.h>
#include <xplorSimulation.h>
#include "prototypes.h"

XplorPot::XplorPot(const String&     instanceName,
			 Simulation* simulation) :
  Pot("XplorPot",instanceName), potIndex(-1)
{
 using CDS::exception;

 if ( simulation->type() != String("XplorSimulation") )
   throw exception("XplorPot::XplorPot: requires an XplorSimulation.");

 sim = (XplorSimulation*)simulation;

 String name = instanceName;
 name.upcase();
 instanceName_ = name;
 if ( name=="ALL" || name=="XPLOR" )
   instanceName_ = "XPLOR";
 else {
   for (int i=0 ; i<sim->xplorVars()->nenert ; i++) {
     String potName(sim->xplorVars()->aner+i*4,4);
     while ( potName.length() &&
	     potName[potName.length()-1] == ' ' ) // remove trailing spaces
       potName.resize( potName.length()-1 );
     
     if ( name == potName ) {
       potIndex = i;
       break;
     }
   }
   if ( potIndex<0 ) {
     String msg = "XplorPot::XplorPot: no such XPLOR potential term: " +
		  name + "\n";
     msg += "valid names are:";
     for (int i=4 ; i<sim->xplorVars()->nenert ; i++) {
       String name = String(sim->xplorVars()->aner+i*4,4);
       if ( name!="    " )
	 msg += " " + name;
     }
     throw exception(msg);
   }
 }

} /* constructor */

EnergyReport 
XplorPot::calcEnergy() 
{
 sim->syncTo();

 int numTerms = sim->xplorVars()->nenert;
 auto_arr<long> qener_save( new long[numTerms] );

 for (int i=0 ; i<numTerms ; i++)
   qener_save[i] = sim->xplorVars()->qener[i];
   
 
 if ( potIndex<0 ) 
   sim->xplorVars()->qener[sim->scriptingIndex()] = 0;
 else {
   for (int i=0 ; i<numTerms ; i++)
     sim->xplorVars()->qener[i] = 0;
   sim->xplorVars()->qener[potIndex] = 1;
 }

 FORTRAN(energy)(sim->xplorVars()->dx,
		 sim->xplorVars()->dy,
		 sim->xplorVars()->dz);

 float_type calcedEnergy=0.;
 for (int i=0 ; i<sim->xplorVars()->nenert ; i++) {
   String name = String(sim->xplorVars()->aner+i*4,4);
   if ( sim->xplorVars()->qener[i] && name!="    " )
     calcedEnergy += sim->xplorVars()->renr[i];
 }

 for (int i=0 ; i<sim->xplorVars()->nenert ; i++)
   sim->xplorVars()->qener[i] = qener_save[i];

 return EnergyReport(instanceName(),calcedEnergy,0.);

} /* calcEnergy */

  //
  // calculate energy and derivatives wrt atomic positions. 
  // Adds these derivs to the values in the DerivList that's passed in.
  // 

EnergyReport 
XplorPot::calcEnergyAndDerivs(DerivList& derivs)
{
 EnergyReport ret = calcEnergy();

 DerivList::VectorVec3& deriv = derivs[ sim ];

 for (int i=0 ; i < sim->xplorVars()->natom ; i++)
   deriv(i) += Vec3(sim->xplorVars()->dx[i],
		    sim->xplorVars()->dy[i],
		    sim->xplorVars()->dz[i]);

 return ret;
} /* calcEnergyAndDerivs */

#endif
