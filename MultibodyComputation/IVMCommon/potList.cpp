
#include <potList.h>

#include <cdsExcept.h>
#include <derivList.h>
#include <simulationWorld.h>
#include <cdsIomanip.h>


PotList::PotList(const String name)
  : Pot("PotList",name)
{
}

static int 
byNameComparer(const rc_Pot& a,
	       const rc_Pot& b) 
{
 int ret = 1;
 if ( String(a->instanceName()) < b->instanceName() )
   ret = -1;
 else if ( String(a->instanceName()) == b->instanceName() )
   ret = 0;

 return ret;
} /* ByNameComparer */


void 
PotList::add(rc_Pot& pot) 
{ 
 if ( pot.ptr() == this )
   throw CDS::exception("PotList::add: can't add self to list.");

 if ( String(pot->instanceName()).length() == 0 )
   throw CDS::exception("PotList::add: term has null name.");

 //
 // make sure we can't add the same potential twice
 //
 bool present=0;
 for (int i=0 ; i<list_.size() ; i++)
   if ( list_[i]->instanceName() == String(pot->instanceName()) ) {
     present=1;
     break;
   }

 if ( !present ) {
   list_.append(pot);
   list_.sort( byNameComparer );
 }
} /* add */

void 
PotList::remove(const String& instanceName)
{

 for (int i=0 ; i<list_.size() ; i++)
   if ( instanceName == list_[i]->instanceName() ) {
     list_.remove(i);
     return;
   }

 throw CDS::out_of_range(String("PotList::remove: no term named ") + 
			 instanceName);

} /* remove */

rc_Pot 
PotList::byName(const String& instanceName)
{

 for (int i=0 ; i<list_.size() ; i++)
   if ( instanceName == list_[i]->instanceName() ) {
     return list_[i];
   }

 throw CDS::out_of_range(String("PotList::byName: no term named ") + 
			 instanceName);

} /* remove */

void 
PotList::removeAll() 
{
 list_.resize(0);
} /* removeAllPots */

EnergyReport
PotList::calcEnergy() 
{
 float_type totEner = 0.0;
 clearEnergyReports();

 for (int i = 0 ; i < list_.size() ; i++) {
   double oldScale = list_[i]->scale();
   list_[i]->setScale( list_[i]->scale() * scale() );
   EnergyReport result = list_[i]->calcEnergy();
   list_[i]->setScale( oldScale );
   addEnergyReport(result);
   totEner += result.energy;
 }

 return EnergyReport(instanceName(),totEner,0.);
} /* calcEnergy */

EnergyReport
PotList::calcEnergyAndDerivs(DerivList& derivList)
{
 float_type totEner = 0.0;
 clearEnergyReports();
 // derivList.clear();

 for (int i = 0 ; i < list_.size() ; i++) {
   double oldScale = list_[i]->scale();
   list_[i]->setScale( list_[i]->scale() * scale() );
   if ( SimulationWorld::world()->logLevel() >=
	SimulationWorld::LOG_DEBUG)
     cout << "PotList::calcEnergyAndDerivs: calling "
	  << list_[i]->instanceName() << endl;
   EnergyReport result = list_[i]->calcEnergyAndDerivs( derivList );
   if ( SimulationWorld::world()->logLevel() >=
	SimulationWorld::LOG_DEBUG)
     cout << "\tresult: " << setprecision(16) << result.energy << endl;   
   list_[i]->setScale( oldScale );
   addEnergyReport(result);
   totEner += result.energy;
 }

 return EnergyReport(instanceName(),totEner,0.);
} /* calcEnergyAndDerivs */

String
PotList::showReport() const
{
 return String("FIX: write me");
} /* showReport */
