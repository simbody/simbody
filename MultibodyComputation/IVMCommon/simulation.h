
#ifndef __simulation_hh__
#define __simulation_hh__

#include <cdsList.h>
#include <cdsMap.h>
#include <cdsString.h>
#include "cdsVec3.h"
#include <atom.h>
#include <cdsPair.h>
#include <cdsVector.h>


//
// Simulation class
//
// An abstract base class for more specific types of simulations.
//
//
// Simulations provide the following:
//
// 1) structure information, including atom information and connectivity.
//
// 2) certain global constants, such that units might be changed.
//
// 3) a method which allows selection of a subset of atoms 
//    (returns a list of atom internal indexes).
//
// 4) a mechanism to calculate and query energy terms associated with atoms
//    within the Simulation.
//
// External to the Simulation: 
//    dynamics and refinement methods
//    certain energy terms
//
// note concerning indexed accessors: all the given index variables are
//  named in such a fashion that the range can be checked by SWIG. Please
//  obey this convention when deriving from this class.
//
//

class PublicIVM;
class ModifiedBase;

class Simulation {

public:
  
  //
  // constructor/destructor
  //

  Simulation(const char* name);

  virtual ~Simulation();

  //
  // accessors
  //

  int id() const { return id_; }

  const char* name() const { return name_; }
  virtual const char* type() const { return "Simulation"; }

  //CDS 11/18/02
  //void addAtom(const Atom& atom) { atomList_.append(atom); }
  //void clearAtomList()           { atomList_.resize(0); }

  void addBondedPair(Pair newVal);
  void clearBondedPairList();

  int numAtoms() const { return numAtoms_; }
  int numBonds() const { return numBonds_; }

  //  const Atom& atomByID(int i)     const { return Atom(this,i); } 
  virtual const Pair& bondPairByID(int bondIndex) const 
  { return bondPairList_[bondIndex]; }

  Atom atomByID(int index) const  { return Atom(this,index); } 
  Atom atomByID(int index)        { return Atom(this,index); } 

  // return string identifying atom by id
  virtual String atomString(int index) const; 



  //
  // routine for initialization of IVM
  //  overridden if necessary
  //
  virtual void registerCallbacks(PublicIVM*) {}

  //
  // access to atom properties
  //
  // Setters are virtual to permit subclasses to do useful
  // things (ala xplorSimulation).  Getters are NOT virtual
  // to permit fast access 
  //

  virtual CDSVector<CDSVec3> atomPosArr() const { return atomPosList_; }
  virtual CDSVector<CDSVec3> atomVelArr() const { return atomVelList_; }
  virtual void setAtomPosArr(const CDSVector<CDSVec3>&);
  virtual void setAtomVelArr(const CDSVector<CDSVec3>&);

  virtual void   setAtomPos(int index, const CDSVec3& newVal);
  virtual void   setAtomVel(int index, const CDSVec3& newVal);
  virtual void  setAtomMass(int index, const float_type& newVal);
  virtual void  setAtomFric(int index, const float_type& newVal);
  virtual void  setAtomCharge(int index, const float_type& newVal);

  virtual void setSegmentName(int index, const char* newVal);
  virtual void setResidueName(int index, const char* newVal);
  virtual void  setResidueNum(int index, const int newVal);
  virtual void    setAtomName(int index, const char* newVal);
  virtual void    setChemType(int index, const char* newVal);
    
  virtual const CDSVec3&        atomPos(int index) const;
  virtual const CDSVec3&        atomVel(int index) const;
  virtual const float_type& atomMass(int index) const;
  virtual const float_type& atomFric(int index) const;
  virtual const float_type& atomCharge(int index) const;
  virtual const char* segmentName(int index) const;
  virtual const char* residueName(int index) const;
  virtual int         residueNum(int index) const;
  virtual const char* atomName(int index) const;
  virtual const char* chemType(int index) const;

  void resizeAtomArrays(const int size);

  //
  // when a Simulation property is updated, all dependent objects are
  // notified of the fact.
  //
  void addDependent(ModifiedBase* dependent) const;
  void removeDependent(ModifiedBase* dependent) const;
  void markAsModified();
  

  //
  // selector--accepts a selection as a char* 
  //

  virtual CDSList< int > select(const char*) const = 0;


  virtual float_type kineticEnergy();

  //
  // thread barrier, etc. call this after modifying coords, etc.
  //
  virtual void sync() {}

  //global default simulation
  static Simulation* currentSimulation() { return currentSimulation_; }
  static void makeCurrent(const Simulation*);

  static int numSimulations();
  static Simulation* simulationByID(int simID);



protected:

  //
  // instance vbls
  //

  int id_;
  String name_;

  int numAtoms_;
  int numBonds_;
  CDSList< String >     atomNameList_;
  CDSList< String >     residueNameList_;
  CDSList< int >        residueNumList_;
  CDSList< String >     segmentNameList_;
  CDSList< String >     chemTypeList_;
  CDSList< float_type > atomMassList_;
  CDSList< float_type > atomFricList_;
  CDSList< float_type > atomRadiusList_;
  CDSList< float_type > atomChargeList_;

  CDSVector< CDSVec3 >       atomPosList_;
  CDSVector< CDSVec3 >       atomVelList_;

  CDSList< Pair >          bondPairList_;
  mutable CDSList< ModifiedBase* > dependentList_;

  static Simulation* currentSimulation_;

};




//
// inline function declarations
//

inline const CDSVec3 &
Simulation::atomPos(int i) const
{
 return atomPosList_[i];
}

inline const CDSVec3 &
Simulation::atomVel(int i) const
{
 return atomVelList_[i] ;
}

inline const float_type &
Simulation::atomMass(int i) const
{
 return atomMassList_[i];
}

inline const float_type &
Simulation::atomFric(int i) const
{
 return atomFricList_[i];
}

inline const char*  
Simulation::segmentName(int i) const
{
 return segmentNameList_[i];
}

inline const char*
Simulation::residueName(int i) const
{
 return residueNameList_[i];
}

inline const char*
Simulation::chemType(int i) const
{
 return chemTypeList_[i];
}

inline int
Simulation::residueNum(int i) const
{
 return residueNumList_[i];
}

inline const char*
Simulation::atomName(int i) const
{
 return atomNameList_[i];
}

inline const float_type&
Simulation::atomCharge(int i) const
{
 return atomChargeList_[i];
}

//inline const float_type &
//Simulation::atomCharge(const int i) const
//{
// return charge();
//}

#endif /* __simulation_hh__ */
