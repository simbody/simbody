#ifdef NOTDEF

#include "xplorSimulation.h"
#include <simulationWorld.h>
#include "xplorWrap.h"
#include <xplorVars.h>
#include <derivList.h>
#include "prototypes.h"
#include <cdsString.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

//
// default initialization of the singleton XplorSimulation instance
//

XplorSimulation* XplorSimulation::simulation_=0;
static int       refCnt;   //inc'ed by initSimulation

//
// constructor
//

XplorSimulation::XplorSimulation() :
  Simulation("XplorSimulation"),
  wrap_( new XplorWrap(this) ),
  modified(0),
  modified_atom( modified_pos | modified_vel | 
		 modified_mass | modified_atomID )
{
 initFrom();
 scriptingIndex_=-1;
 for (int i=0 ; i<xplorVars()->nenert ; i++)
   if ( String(xplorVars()->aner+i*4,4) == "SCRI" )
     scriptingIndex_ = i;
 
 if ( scriptingIndex_<0 )
   throw CDS::exception(String("XplorSimulation::XplorSimulation: ") +
			"could not find XPLOR SCRI energy term");
} /* constructor */

//
// destructor
//

XplorSimulation::~XplorSimulation()
{ 
 delete wrap_;
}

XplorVars* 
XplorSimulation::xplorVars()
{
 return XplorWrap::xplorVars();
} /* xplorVars */

void
XplorSimulation::setAtomPosArr(const CDSVector<CDSVec3>& arr)
{
 modified |= modified_pos;
 Simulation::setAtomPosArr(arr);
} /* setAtomPosArr */

void
XplorSimulation::setAtomVelArr(const CDSVector<CDSVec3>& arr)
{
 modified |= modified_vel;
 Simulation::setAtomVelArr(arr);
} /* setAtomPosArr */



//
// accessors
//

void
XplorSimulation::setAtomPos(      int   i, 
			    const CDSVec3& newVal) 
{
 modified |= modified_pos;
 Simulation::setAtomPos(i, newVal);
}

void
XplorSimulation::setAtomVel(      int   i, 
			    const CDSVec3& newVal)
{
 modified |= modified_vel;
 Simulation::setAtomVel(i, newVal);
}

void
XplorSimulation::setAtomMass(      int         i, 
			     const float_type& newVal)
{
 modified |= modified_mass;
 Simulation::setAtomMass(i, newVal);
}

void
XplorSimulation::setAtomFric(      int         i, 
			     const float_type& newVal) 
{
 modified |= modified_fric;
 Simulation::setAtomFric(i, newVal);
}

void
XplorSimulation::setAtomCharge(      int         i, 
			       const float_type& newVal) 
{
 modified |= modified_charge;
 Simulation::setAtomCharge(i, newVal);
}

void
XplorSimulation::setSegmentName(      int i, 
				const char* newVal) 
{
 modified |= modified_atomID;
 Simulation::setSegmentName(i, newVal);
}

void
XplorSimulation::setResidueName(      int i, 
				const char* newVal) 
{
 modified |= modified_atomID;
 Simulation::setResidueName(i, newVal);
}

void
XplorSimulation::setResidueNum(      int i, 
			       const int newVal) 
{
 modified |= modified_atomID;
 Simulation::setResidueNum(i, newVal);
}

void
XplorSimulation::setAtomName(      int i, 
			     const char* newVal) 
{
 modified |= modified_atomID;
 Simulation::setAtomName(i, newVal);
}

void
XplorSimulation::setChemType(      int i, 
			     const char* newVal) 
{
 modified |= modified_atomID;
 Simulation::setChemType(i, newVal);
}

//void
//XplorSimulation::setFullName(const int i, 
//			     const char* newVal) 
//{
// modified |= modified_atomID;
// Simulation::setFullName(i, newVal);
//}

//
// Makes sure that the number of atoms on the scripting side is
// the same as the number of atoms on the xplor side.  Also
// does the same thing for the bondedPairs.  
// Intended to be used ONLY within sync... calls, since it does
// NOT actually sync the atom data across
//

class XplorRandomSeeder : public RandomSeeder {
  const XplorWrap* wrap;
public:
  XplorRandomSeeder(const XplorWrap* wrap) : wrap(wrap) {}
  virtual ~XplorRandomSeeder() {}
  virtual void setSeed(int seed) { wrap->setRandomSeed(seed); }
};


void
XplorSimulation::initFrom()
{

 SimulationWorld::init(xplorVars()->timeFac  ,
		       xplorVars()->kBoltzman,
		       (SimulationWorld::LogLevel)xplorVars()->wrnlev);

 SimulationWorld::world()->registerRandomSeeder(new XplorRandomSeeder(wrap()));

 //
 // if we don't have the same number of atoms
 // as xplor has, dump our atom list and 
 // create a new one.  Does NOT sync up the
 // data afterward
 //

 if ( xplorVars()->natom != numAtoms() ) {

   resizeAtomArrays(xplorVars()->natom);

   syncFrom();

 }

 //
 // if our bond list doesn't have the same
 // number of entries as xplor's, dump ours
 // and copy in xplor's
 //

 if ( xplorVars()->nbond != numBonds() ) {

   clearBondedPairList();
   
   for (int i=0 ; i < xplorVars()->nbond ; i++)
     //
     // note that atom indices obtained from the
     // xplor ib and jb arrays are switched to zero-offset.
     //
     addBondedPair( Pair(xplorVars()->ib[i]-1,xplorVars()->jb[i]-1) );
 }

}  /* initFrom */



//
// the proper top-level interface 
//

void
XplorSimulation::syncFrom()
{
 syncAtomIDFrom();
 syncPosFrom();
 syncVelFrom();
 syncFricFrom();
 syncChargeFrom();
 syncMassFrom();
 sync();  
}

static String
fromXplorString(const char* start,
		      int   len  )
  // create String from char* with max length len
  // remove trailing space
{
 while ( len>0 && start[len-1]==' ' ) len--;
 return String(start,len);
} /* fromXplorString */

void
XplorSimulation::syncAtomIDFrom() 
{
 initFrom();
  
 //
 // the xplor string arrays don't end with \0s, so I have
 // to break them up myself
 //
  
 for (int i=0; i < xplorVars()->natom; i++) {
   segmentNameList_[i] = fromXplorString(xplorVars()->segName+i*4,4);
   residueNameList_[i] = fromXplorString(xplorVars()->resName+i*4,4);
   residueNumList_[i]  = atoi( String(xplorVars()->resNum+i*4,4) );
   atomNameList_[i]    = fromXplorString(xplorVars()->atomName+i*4,4);
   chemTypeList_[i]    = fromXplorString(xplorVars()->chemType+i*4,4);
 }
}

void
XplorSimulation::syncFricFrom() 
{
 initFrom();

 for (int i=0; i < xplorVars()->natom; i++)
   atomFricList_[i] = xplorVars()->fbeta[i];
  
}

void
XplorSimulation::syncChargeFrom() 
{
 initFrom();

 for (int i=0; i < xplorVars()->natom; i++)
   atomChargeList_[i] = xplorVars()->charge[i];
}

void
XplorSimulation::syncMassFrom() 
{
 initFrom();
  
 for (int i=0; i < xplorVars()->natom; i++)
   atomMassList_[i] = xplorVars()->mass[i];
  
}

void
XplorSimulation::syncPosFrom()
{
 initFrom();

 for (int i=0 ; i < xplorVars()->natom ; i++)
   if ( xplorVars()->x[i]<9998.0 )
     atomPosList_[i] =  CDSVec3(xplorVars()->x[i],
			    xplorVars()->y[i],
			    xplorVars()->z[i]);
   else
     atomPosList_[i] = CDSVec3(Atom::INVALID_COORD,
			   Atom::INVALID_COORD,
			   Atom::INVALID_COORD);
} /* syncPosFrom */

void
XplorSimulation::syncVelFrom()
{
 initFrom();

 for (int i=0 ; i < xplorVars()->natom ; i++)
   atomVelList_[i] = CDSVec3(xplorVars()->xv[i],
			 xplorVars()->yv[i],
			 xplorVars()->zv[i]);

}

//
// convenience functions for use in syncTo, below
//

//
// pad a string with trailing spaces to length 4
//

static const char*
pad4(const char* s)
{
 static char* ret = new char[5];
 sprintf(ret, "%-4s", s);
 return ret;
}

//
// convert an int to a padded string with trailing spaces to length 4
//

static const char*
pad4(int i)
{
 static char* ret = new char[5];
 sprintf(ret, "%-4d", i);
 return ret;
} /* pad4 */


//
// copy XplorSimulation's values to xplor
//

void
XplorSimulation::syncTo()
{
 assert( xplorVars()->natom==atomNameList_.size() );


 if ( modified & modified_atomID )
   for (int i=0; i < xplorVars()->natom; i++) {

     //
     // xplorVars are just character arrays, with no trailing \0s,
     // so I have to use memcpy to write to them.  I have to use
     // the sprintf call to ensure that names that are shorter than
     // 4 chars are padded with spaces instead of random junk.
     // Note that the residue num array is actually character-typed, 
     // so I need a sprintf there as well.
     //

       
     memcpy(xplorVars()->segName+i*4,  pad4(segmentNameList_[i]), 4);
     memcpy(xplorVars()->resName+i*4,  pad4(residueNameList_[i]), 4);
     memcpy(xplorVars()->resNum+i*4,   pad4(residueNumList_[i] ), 4);
     memcpy(xplorVars()->atomName+i*4, pad4(atomNameList_[i]   ), 4);
     memcpy(xplorVars()->chemType+i*4, pad4(chemTypeList_[i]   ), 4);
   }

 if ( modified & modified_pos )
   for (int i=0 ; i < xplorVars()->natom ; i++) {
     xplorVars()->x[i] = atomPosList_[i].x();
     xplorVars()->y[i] = atomPosList_[i].y();
     xplorVars()->z[i] = atomPosList_[i].z();
   }

 if ( modified & modified_vel )
   for (int i=0 ; i < xplorVars()->natom ; i++) {
     xplorVars()->xv[i] = atomVelList_[i].x();
     xplorVars()->yv[i] = atomVelList_[i].y();
     xplorVars()->zv[i] = atomVelList_[i].z();
   }

  
 if ( modified & modified_mass )
   for (int i=0 ; i < xplorVars()->natom ; i++)
     xplorVars()->mass[i] = atomMassList_[i];

 if ( modified & modified_fric )
   for (int i=0 ; i<xplorVars()->natom ; i++)
     xplorVars()->fbeta[i] = atomFricList_[i];
  
 if ( modified & modified_charge )
   for (int i=0 ; i<xplorVars()->natom ; i++)
     xplorVars()->charge[i] = atomChargeList_[i];

 modified = 0;
} /* syncTo */

//
// Atom selector, which just calls xplor's version
//

CDSList< int > 
XplorSimulation::select(const char* sel) const
{ 
 return wrap_->select(sel); 
} /* select */ 


//
// calcNonXplorEnergy is called ONLY by xplor's
// FORTRAN-side main energy routine. 
//
// If xplor is in control, the scripting side needs to 
// be synched up from xplor before it can calculate
// the scripting energies.
//
// We assume that ONLY the positions change!
//
// If xplor is in control, the first time this is 
// called, the number of atoms on the script side 
// will be zero.  In syncPosFrom, this will fail
// the test about the number of atoms, forcing a 
// complete sync the first time.  
//
// But if the scripting side is in control, 
// then this is a bit wasteful, since there's no 
// need to do a syncPosFrom (since the coords
// were just copied to xplor by the syncTo at
// the top of calcEnergyAndDeriv()).  
//

float_type
XplorSimulation::calcNonXplorEnergyAndDeriv(DerivList& derivList) 
{
 //called before other xplor energy terms so that we don't need to 
 // do a syncDerivFrom
 syncPosFrom();
 // float_type scriptingEner = Simulation::calcEnergyAndDerivs(derivList);
 float_type scriptingEner = potList_.calcEnergyAndDerivs(derivList).energy;
 syncTo();

 return scriptingEner;
}


//
// The energy routine that is actually called by XPLOR
//

extern "C" double
FORTRAN(scripting_energy)()
{
 //
 // make sure we don't calculate energy if the singleton
 // instance of XplorSimulation hasn't been initialized
 //
 XplorSimulation* sim = XplorSimulation::simulation();
 if ( sim == 0 )
   return 0.0;
 if ( !sim->xplorVars()->qener[sim->scriptingIndex()] )
   return 0.;

 static DerivList derivList;
 derivList.init( sim );
 
 float_type nonXplorEner = 
   sim->calcNonXplorEnergyAndDeriv(derivList);

 const CDSVector<CDSVec3>& deriv = derivList[sim];
 XplorVars* xplorVars = sim->xplorVars();
 for (int i=0 ; i < xplorVars->natom ; i++) {
   xplorVars->dx[i] += deriv(i).x();
   xplorVars->dy[i] += deriv(i).y();
   xplorVars->dz[i] += deriv(i).z();
 }

 return nonXplorEner;
} /* scripting_energy */

//
// initialization of the SINGLETON instance of XplorSimulation
//

void
XplorSimulation::initSimulation()
{
 if ( simulation_ ) 
   return;

 simulation_ = new XplorSimulation();
 refCnt++;

 currentSimulation_ = simulation_; //global default simulation
} /* initSimulation */

void
XplorSimulation::deleteSimulation()
{
 if ( !simulation_ )
   return;

 refCnt--;
 if ( refCnt==0 )
   delete simulation_;
 simulation_=0;
} /* deleteSimulation */

#endif
