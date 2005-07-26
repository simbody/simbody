#ifdef NOTDEF


#ifndef __scriptingsimulation_hh__
#define __scriptingsimulation_hh__

#include <simTypes.h>
#include <simulation.h>

//
// ScriptingSimulation
//
// This class summarizes the Simulation class heirarchy's
// interface that should be accessible to the scripting level.
// 
// Mostly, this consists of accessor functions.  These accessors
// include more expensive validity tests than the low-level 
// programmatic interfaces in the Simulation class heirarchy.
//

class ScriptingSimulation {

  //inaccessible
  ScriptingSimulation();
  ScriptingSimulation(const ScriptingSimulation&);
  ScriptingSimulation& operator=(const ScriptingSimulation);

public:

  //
  // constructor
  //

  ScriptingSimulation(Simulation* aSimulationPtr);

//  // conversion to Simulation*
//  operator Simulation*() {return sim;}
  Simulation* simulation() { return sim; }

  //
  // accessors
  //

  const char* name() const;

  int numAtoms() const;
  int numBonds() const;

  //const Atom& atomByID(int i);
  Pair bondPairByID(int i);

  CDSList<Vec3> atomPosArr() const             { return sim->atomPosArr(); }
  void setAtomPosArr(const CDSList<Vec3>& arr) { sim->setAtomPosArr(arr); }

  void   setAtomPos(const int i, const Vec3& newVal);
  void   setAtomVel(const int i, const Vec3& newVal);
  void  setAtomMass(const int i, const float_type newVal);
  void  setAtomFric(const int i, const float_type newVal);
  void  setAtomCharge(const int i, const float_type newVal);
  
  void setSegmentName(const int i, const char* newVal);
  void setResidueName(const int i, const char* newVal);
  void  setResidueNum(const int i, const int newVal);
  void    setAtomName(const int i, const char* newVal);
  //  void    setFullName(const int i, const char* newVal);
  
  Vec3        atomPos(const int i) const;
  Vec3        atomVel(const int i) const;
  float_type  atomMass(const int i) const;
  float_type  atomFric(const int i) const;
  float_type  atomCharge(const int i) const;
  
  const char* segmentName(const int i) const;
  const char* residueName(const int i) const;
  int         residueNum(const int i) const;
  const char* atomName(const int i) const;
  //  const char* fullName(const int i) const { return sim->fullName(i); }

 
  //
  // selector
  //

  CDSList< int > select(const char*);


  static Simulation* currentSimulation();

private:

  //
  // tests
  //

  void checkIndexRange(const int i) const;

  //
  // instance vbls
  //

  Simulation* sim;

};


#endif /* __scriptingsimulation_hh__ */

#endif
