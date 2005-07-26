
#ifndef __xplorsimulation_hh__
#define __xplorsimulation_hh__

#include "simulation.h"
#include <potList.h>
#include <cdsVector.h>

class XplorVars;
class XplorWrap;

//
// XplorSimulation
//
// A subclass of Simulation that is tied to xplor.
// 
// It maintains its own copy of the most important xplor variables
// (position, velocity arrays, etc)
//
// It also maintains its own copy of the xplor wrapper functions
// to be able to make calls to select, among other things
//
// Can sync its copies of the xplor variables back and forth 
// with xplor
//
// Note that the GLOBAL xplorSimulation and xplorVars ptrs
// are defined at the bottom of this file!
//


class XplorSimulation : public Simulation {

public:
  
  // 
  // constructor/destructor
  //

  XplorSimulation();
  virtual ~XplorSimulation();

  virtual const char* type() const { return "XplorSimulation"; }

  //
  // these accessor functions are overridden to ensure that
  // the proper modified... flags are set
  //

  virtual void setAtomPosArr(const CDSVector<Vec3>&);
  virtual void setAtomVelArr(const CDSVector<Vec3>&);

  virtual void   setAtomPos(      int   i, 
			    const Vec3& newVal);
  virtual void   setAtomVel(      int i, 
			    const Vec3& newVal);
  virtual void  setAtomMass(      int i, 
			    const float_type& newVal);
  virtual void  setAtomFric(      int i, 
			    const float_type& newVal);
  virtual void  setAtomCharge(      int i, 
			      const float_type& newVal);

  virtual void setSegmentName(      int i, 
			      const char* newVal);
  virtual void setResidueName(      int i, 
			      const char* newVal);
  virtual void  setResidueNum(      int i, 
			      const int newVal);
  virtual void    setAtomName(      int i, 
			      const char* newVal);
  virtual void    setChemType(      int i, 
			      const char* newVal);

  //
  // returns a const pointer to the xplor wrap
  //

  const XplorWrap* wrap() const { return wrap_; }

  //
  // atom selector
  //

  CDSList< int > select(const char* sel) const;

  
  PotList& potList() { return potList_; }

  XplorVars* xplorVars();

  //
  // Routines to synchronize simulation values to/from xplor
  //

  void initFrom();       //                                             
  void syncFrom(); 	 // initialize/copy instance variables from xplor
  void syncAtomIDFrom(); //                                           
  void syncPosFrom();
  void syncVelFrom();
  void syncFricFrom();
  void syncChargeFrom();
  void syncMassFrom();

  void syncTo();      // copy modified instance variables to xplor

  // index into XPLOR renr array corresponding to scripting term
  int scriptingIndex() { return scriptingIndex_; }

private:

  // 
  // instance vbls
  //

  XplorWrap* wrap_;
  unsigned modified;
  unsigned modified_atom;
  int scriptingIndex_;
  PotList potList_;

  //
  // convenient flags for the sync process
  //

  enum { modified_pos    = 1, 
	 modified_vel    = 2, 
	 modified_mass   = 8, 
	 modified_fric   = 16,
	 modified_charge = 32,
	 modified_atomID = 64 };

  //
  // singleton XplorSimulation ptr
  //
  
  static XplorSimulation* simulation_;
public:
  static void initSimulation();
  static void deleteSimulation();
  static XplorSimulation* simulation() { return simulation_; }

  //FIX: should be in private interface
  float_type calcNonXplorEnergyAndDeriv(DerivList&); 
  //FIX: scripting_energy should be made a friend - is this possible?

};


#endif /* __xplorsimulation_hh__ */
