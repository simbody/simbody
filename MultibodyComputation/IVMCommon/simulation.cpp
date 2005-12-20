
#include "simulation.h"

#include <cdsExcept.h>
#include <cdsSStream.h>
#include <cdsIomanip.h>
#include <modified.h>

Simulation* Simulation::currentSimulation_ = 0;
static CDSList<Simulation*>* simulationList=0;
static int lastSimulation=0;

Simulation::Simulation(const char* name) :
  name_(name), numAtoms_(0), numBonds_(0)
{ 
 id_ = lastSimulation;
 lastSimulation++;

 if (!simulationList)
   simulationList = new CDSList<Simulation*>;

 simulationList->append( this );
} /* constructor */

Simulation::~Simulation()
{
 int index = simulationList->getIndex( this );
 if ( index>=0 )
   simulationList->remove(index);
 for (int i=0 ; i<dependentList_.size() ; i++)
   dependentList_[i]->unRegister(this);
} /* destructor */


//
// accessors
//

void 
Simulation::addBondedPair(Pair newVal)
{
 numBonds_++;
 bondPairList_.append(newVal);
 markAsModified();
} /* addBondedPair */

void
Simulation::clearBondedPairList() 
{
 numBonds_=0;
 bondPairList_.resize(0);
 markAsModified();
} /* clearBondedPairList */



//
// convenience accessors
//


void
Simulation::setAtomPosArr(const CDSVector<CDSVec3>& arr)
{
 if (arr.size() != numAtoms())
   throw CDS::exception("FIX: setAtomPosArr: array should have size numAtoms");

 atomPosList_ = arr;
 markAsModified();
} /* setAtomPosArr */

void
Simulation::setAtomVelArr(const CDSVector<CDSVec3>& arr)
{
 if (arr.size() != numAtoms())
   throw CDS::exception("FIX: setAtomVelArr: array should have size numAtoms");

 atomVelList_ = arr;
 markAsModified();
} /* setAtomVelArr */


void 
Simulation::setAtomPos(      int i, 
		       const CDSVec3& newVal)
{
 atomPosList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setAtomVel(      int i, 
		       const CDSVec3& newVal)
{
 atomVelList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setAtomMass(      int i, 
			const float_type& newVal) 
{
 atomMassList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setAtomFric(      int i, 
			const float_type& newVal) 
{
 atomFricList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setAtomCharge(      int i, 
			  const float_type& newVal)
{
 atomChargeList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setSegmentName(      int i, 
			   const char* newVal) 
{
 segmentNameList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setResidueName(      int i, 
			   const char* newVal) 
{
 residueNameList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setChemType(      int i, 
			const char* newVal) 
{
 chemTypeList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setResidueNum(      int i, 
			  const int newVal) 
{
 residueNumList_[i] = newVal;
 markAsModified();
}

void 
Simulation::setAtomName(      int i, 
			const char* newVal) 
{
 atomNameList_[i] = newVal;
 markAsModified();
}

String
Simulation::atomString(int id) const
{
 StringStream ret;
 ret << id << ' ';
 if (id<0 || id>=numAtoms())
   ret << "(unknown)";
 else
   ret << setw(4) << residueNum(id)  << ' ' 
       << setw(4) << residueName(id) << ' ' 
       << setw(4) << atomName(id);
 ret << ends;
 return ret.str();
} /* atomString */


float_type 
Simulation::kineticEnergy()
{
 float_type ret =0.;
 for (int i=0 ; i<numAtoms() ; i++)
   ret += atomMass(i) * dot( atomVel(i),atomVel(i) );
 return 0.5*ret;
} /* kineticEnergy */

void
Simulation::makeCurrent(const Simulation* sim)
{
 currentSimulation_ = (Simulation*)sim;
}

int
Simulation::numSimulations()
{
 if (!simulationList) return 0;
 return simulationList->size();
}

Simulation*
Simulation::simulationByID(int id)
{
 if (!simulationList) return 0;
 return (*simulationList)[ id ];
}

void
Simulation::resizeAtomArrays(const int size)
{
 int oldsize = atomNameList_.size();
 if ( size==oldsize ) return;

 markAsModified();

 atomNameList_.resize(size);
 residueNameList_.resize(size);
 residueNumList_.resize(size);
 segmentNameList_.resize(size);
 chemTypeList_.resize(size);
 atomPosList_.resize(size);
 atomVelList_.resize(size);
 atomMassList_.resize(size);
 atomFricList_.resize(size);
 atomRadiusList_.resize(size);
 atomChargeList_.resize(size);

 numAtoms_ = size;

 //
 // atom property innitialization 
 //
 for (int i=oldsize ; i<size ; i++) {
   atomNameList_[i] = String( "none",4,0 );
   residueNameList_[i] = String( "none",4,0 );
   residueNumList_[i]  = -1;
   segmentNameList_[i] = String( "none",4,0 );
   atomPosList_[i] = CDSVec3(Atom::INVALID_COORD, 
			  Atom::INVALID_COORD, 
			  Atom::INVALID_COORD);
   atomVelList_[i] = CDSVec3(0, 0, 0);
   atomMassList_[i] = 0;
   atomFricList_[i] = 0;
   atomRadiusList_[i] = 0;
 }
} /* resizeAtomArrays */

void
Simulation::addDependent(ModifiedBase* dependent) const
{
 dependentList_.append(dependent);
} /* addDependent */

void
Simulation::removeDependent(ModifiedBase* dependent) const
{
 if ( dependentList_.contains(dependent) )
   dependentList_.remove( dependentList_.getIndex(dependent) );
} /* removeDependent */

void
Simulation::markAsModified()
{
 for (int i=0 ; i<dependentList_.size() ; i++)
   dependentList_[i]->modified.set();
} /* markAsModified */
