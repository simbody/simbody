
#ifndef __derivList_hh__
#define __derivList_hh__

#include <cdsList.h>
#include <cdsVector.h>
#include "cdsVec3.h"
#include <mmapAlloc.h>
#include <simulation.h>
#include <atom.h>

class Atom;
class Simulation;

class DerivList {
public:
  typedef CDSVector<CDSVec3,0> VectorVec3;

  CDSList< VectorVec3 > derivList;
  
  DerivList();
  ~DerivList() {}

  //allocate and zero deriv array with given id
  void init(const Simulation*); 

  // zero out all existing arrays
  void clear(); 

  // indexing operators for access to arrays

  VectorVec3& operator[](const Simulation* sim)
  {
   if ( derivList.size() <= sim->id() ||
	sim->numAtoms() != derivList[sim->id()].size() )
     init( sim );
   return derivList[sim->id()];
  }

  const VectorVec3& operator[](const Simulation* sim) const
  {
   return derivList[sim->id()];
  }

//  CDSVec3& operator[](const Atom& anAtom)
//  {
//   return (*this)[anAtom.simulation()][anAtom.index()];
//  }
  CDSVec3& operator[](const Atom& anAtom)
  {
   return (*this)[anAtom.simulation()][anAtom.index()];
  }

  const CDSVec3& operator[](const Atom& anAtom) const
  {
   return (*this)[anAtom.simulation()][anAtom.index()];
  }

};

#endif /* __derivList_hh__ */
