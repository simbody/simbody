
#ifdef NOTDEF

#include "scriptingSimulation.h"
#include <simulation.h>
#include <pot.h>

//
// constructor
//

ScriptingSimulation::ScriptingSimulation(Simulation *aSimulationPtr) {

  sim = aSimulationPtr;
}

//
// tests
//

void 
ScriptingSimulation::checkIndexRange(const int i) const {

  if ( i<0 || i>= sim->numAtoms() )
    throw CDS_NAMESPACE(out_of_range)();
}


//
// accessors
//

const char* 
ScriptingSimulation::name() const { 
  
  return sim->name(); 
}


int 
ScriptingSimulation::numAtoms() const { 

  return sim->numAtoms();
}

int 
ScriptingSimulation::numBonds() const { 

  return sim->numBonds();
}

//const Atom&
//ScriptingSimulation::atomByID(int i)
//{ return sim->atomByID(i); }

Pair
ScriptingSimulation::bondPairByID(int i)
{ return sim->bondPairByID(i); }

void
ScriptingSimulation::setAtomPos(const int i, const CDSVec3& newVal) {

  checkIndexRange(i);
  sim->setAtomPos(i, newVal);
}

void
ScriptingSimulation::setAtomVel(const int i, const CDSVec3& newVal) {

  checkIndexRange(i);
  sim->setAtomVel(i, newVal);
}

void
ScriptingSimulation::setAtomMass(const int i, const float_type newVal) {

  checkIndexRange(i);
  sim->setAtomMass(i, newVal);
}

void
ScriptingSimulation::setAtomFric(const int i, const float_type newVal) {

  checkIndexRange(i);
  sim->setAtomFric(i, newVal);
}

void
ScriptingSimulation::setAtomCharge(const int i, const float_type newVal) {

  checkIndexRange(i);
  sim->setAtomCharge(i, newVal);
}

void 
ScriptingSimulation::setSegmentName(const int i, const char* newVal) {

  checkIndexRange(i);
  sim->setSegmentName(i, newVal);
}

void 
ScriptingSimulation::setResidueName(const int i, const char* newVal) {

  checkIndexRange(i);
  sim->setResidueName(i, newVal);
}

void 
ScriptingSimulation::setResidueNum(const int i, const int newVal) {

  checkIndexRange(i);
  sim->setResidueNum(i, newVal);
}

void 
ScriptingSimulation::setAtomName(const int i, const char* newVal) {

  checkIndexRange(i);
  sim->setAtomName(i, newVal);
}

//void 
//ScriptingSimulation::setFullName(const int i, const char* newVal) {
//
//  checkIndexRange(i);
//  sim->setFullName(i, newVal);
//}


CDSVec3         
ScriptingSimulation::atomPos(const int i) const {

  checkIndexRange(i);
  return sim->atomPos(i);
}

CDSVec3         
ScriptingSimulation::atomVel(const int i) const {

  checkIndexRange(i);
  return sim->atomVel(i);
}

float_type         
ScriptingSimulation::atomMass(const int i) const {

  checkIndexRange(i);
  return sim->atomMass(i);
}

float_type         
ScriptingSimulation::atomFric(const int i) const {

  checkIndexRange(i);
  return sim->atomFric(i);
}

float_type         
ScriptingSimulation::atomCharge(const int i) const {

  checkIndexRange(i);
  return sim->atomCharge(i);
}

const char* 
ScriptingSimulation::segmentName(const int i) const {

  checkIndexRange(i);
  return sim->segmentName(i);
}

const char* 
ScriptingSimulation::residueName(const int i) const {

  checkIndexRange(i);
  return sim->residueName(i);
}

int
ScriptingSimulation::residueNum(const int i) const {

  checkIndexRange(i);
  return sim->residueNum(i);
}

const char* 
ScriptingSimulation::atomName(const int i) const {

  checkIndexRange(i);
  return sim->atomName(i);
}


//
// selector
//

CDSList< int >
ScriptingSimulation::select(const char* sel) {

  return sim->select(sel);
}


Simulation*
ScriptingSimulation::currentSimulation() 
{
 return Simulation::currentSimulation();
}

#endif
