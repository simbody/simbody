#ifdef NOTDEF


#ifndef __ensembleSimulation_hh__
#define __ensembleSimulation_hh__

#include "simulation.h"
#include "atomSel.h"
#include "mmapAlloc.h"
#include "atomSelAction.h"

#include <semaphore.h>


#include <cdsAuto_ptr.h>


//
// EnsembleSimulation
// contains an ensemble of Simulations having the same number of atoms, 
//  topologies, etc.
//
// EnsembleMemberSimulation
//   allows one to access members of the ensembles as simulations.
//
//

class EnsembleMemberSimulation;
class DistributedProduct;
class AllocData;

class EnsembleSimulation : public Simulation {

public:
  
  
  // create size number of subprocesses. 
  // For each one:
  //  1) create semaphore, shared memory region
  //     data to share:
  //      pos, vel, ?
  //  2) then fork
  //  3) each child needs an id to identify itself.
  //  4) master then runs NUM_THREADS number of processes at a time.
  //     the others wait.
  //  5) processes run until need to communicate or share data:
  //      avePos/ aveVel
  //     In these routines, the children block and yield to the master process.
  //  6) output of children is redirected to /dev/null, or a user-determinable
  //      location.
  //  7) Number of atoms cannot change after construction.
  //  8) only one ensemble with size>1 may exist at a time
  //
  EnsembleSimulation(const char*       name,
			   int         size,
			   int         numThreads=-1,
			   Simulation* defaultSimulation=0);

  
  //
  // reaps child processes, clean up shared memory.
  //
  virtual ~EnsembleSimulation();
  
  virtual const char* type() const { return "EnsembleSimulation"; }

  //
  // Synchronize all threads. id is any nonzero integer. It's used to catch
  // barrier mismatches.
  //
  void barrier(const char* module="",
		     int   id=1); 

  //
  // explicitly shutdown threads- not usually needed/desired.
  //
  void shutdown();

  //
  //
  //
  int numThreads() { return numThreads_; }

  //
  // size of ensemble.
  //
  int size() const {return size_;}
  //
  // simulation upon which the ensemble is based.
  //
  Simulation* subSim() const { return subSim_; }
  //
  // return Simulation corresponding to ensemble member.
  //
  EnsembleMemberSimulation* members(int memberIndex)
  { return members_[memberIndex].get(); }
  const EnsembleMemberSimulation* members(int memberIndex) const
  { return members_[memberIndex].get(); }
  //
  // member for this process
  //
  EnsembleMemberSimulation* member() { return member_; }
  const EnsembleMemberSimulation* member() const { return member_; }
  //
  // weight of member i in the ensemble
  //
  float_type weight(int memberIndex) const { return weights_[memberIndex]; }

  virtual float_type kineticEnergy();

  virtual void registerCallbacks(PublicIVM*);

  //
  // special definition of vector dot product distributed across ensemble.
  //
  template<class VEC>
  float_type vecProd(const VEC& v1, const VEC& v2);


  // need: accessor for avePos, aveVel, 
  
  
  // these accessors are forwarded to the memberSimulation
  //
  // accessors
  //
  virtual const Pair& bondPairByID(int bondIndex) const;

  // return string identifying atom by id
  virtual CDSString atomString(int index) const;

  //
  // custom selector
  //
  CDSList< int > select(const char* str) const;

  //
  // access to atom properties
  //

  CDSVector<CDSVec3> meanAtomPosArr() const;
  virtual CDSVector<CDSVec3> atomPosArr() const;
  virtual CDSVector<CDSVec3> atomVelArr() const;
  virtual void setAtomPosArr(const CDSVector<CDSVec3>&);
  virtual void setAtomVelArr(const CDSVector<CDSVec3>&);

  virtual void setAtomPos(int index, const CDSVec3& newVal);
  virtual void setAtomVel(int index, const CDSVec3& newVal);
  virtual void setAtomMass(int index, const float_type& newVal);
  virtual void setAtomFric(int index, const float_type& newVal);
  virtual void setAtomCharge(int index, const float_type& newVal);
  virtual void setSegmentName(int index, const char* newVal);
  virtual void setResidueName(int index, const char* newVal);
  virtual void setResidueNum(int index, const int newVal);
  virtual void setAtomName(int index, const char* newVal);
  virtual void setChemType(int index, const char* newVal);
    
  virtual const CDSVec3& atomPos(int index) const;
  virtual const CDSVec3& atomVel(int index) const;
  virtual const float_type& atomMass(int index) const;
  virtual const float_type& atomFric(int index) const;
  virtual const float_type& atomCharge(int index) const;
  virtual const char* segmentName(int index) const;
  virtual const char* residueName(int index) const;
  virtual int         residueNum(int index) const;
  virtual const char* atomName(int index) const;
  virtual const char* chemType(int index) const;

  virtual void sync();
  void resize();   //reset size from subSim


  //
  // memory allocation routines: to get shared memory after a fork.
  // All members must call this function.
  //
  struct SharedAlloc {
    static void* alloc(size_t);
    static void free(void*);
  };

  static EnsembleSimulation* currentSimulation() {
   return currentSimulation_;
  }

  // return a size-one EnsembleSimulation based on sim
  static EnsembleSimulation* sizeOneSimulation(Simulation* sim=0);



protected:

  //
  // instance vbls
  //

  static EnsembleSimulation* currentSimulation_;

  Simulation* subSim_;
  int size_;
  int numThreads_;
  CDSList<float_type> weights_;
  mutable CDSList< CDS::auto_ptr<EnsembleMemberSimulation> > members_;
  EnsembleMemberSimulation* member_;

  CDSList< CDS::auto_ptr<CDSVector<float_type,0,SharedAlloc> > > sharedReports;

  bool deleted;

public:

#ifndef SWIG
  CDS::auto_ptr<DistributedProduct> distributedProduct;
  CDS::auto_ptr<AllocData>          allocData;
#endif

  friend class EnsembleMemberSimulation;
};


//
//
//

class EnsembleMemberSimulation : public Simulation {
  // constructed by EnsembleSimulation
  EnsembleMemberSimulation(const char*               name,
				 EnsembleSimulation* sim,
			   const int                 memberIndex);
public:
  typedef CDSVector< CDSVec3,0,EnsembleSimulation::SharedAlloc > VectorVec3;
  
  //
  // constructor/destructor
  //

  //EnsembleMemberSimulation();
  //EnsembleMemberSimulation(const char*               n,
  //				 EnsembleSimulation* s,
  //			   const int                 i);

  virtual ~EnsembleMemberSimulation();

  //
  //
  pid_t pid() { return memberData_->pid; }

  //
  //
  int barrierID()    { return memberData_->barrierID; }
  int barrierCnt()   { return memberData_->barrierCnt; }
  void barrierIncr() { memberData_->barrierCnt++; }

  //
  //
  //
  void sleep();
  void wake();
  //bool isAwake();

  //
  // 
  //
  int memberIndex() const { return memberIndex_; }
  //
  // weight for this member
  //
  float_type weight() const { return ensembleSim()->weight(memberIndex_); }
  //
  // 
  //
  EnsembleSimulation* ensembleSim() const { return ensembleSim_; }
  //
  // 
  //
  Simulation* subSim() const { return ensembleSim()->subSim(); }

  //resync atom/bond info from subSim
  void resize();

  const Pair& bondPairByID(int bondIndex) const 
  { return subSim()->bondPairByID(bondIndex); }

  virtual CDSString atomString(int index) const 
  { return subSim()->atomString(index); }


  virtual CDSVector<CDSVec3> atomPosArr() const;
  virtual CDSVector<CDSVec3> atomVelArr() const;
  virtual void setAtomPosArr(const CDSVector<CDSVec3>&);
  virtual void setAtomVelArr(const CDSVector<CDSVec3>&);

  virtual void setAtomPos(int index, const CDSVec3& newVal);
  virtual void setAtomVel(int index, const CDSVec3& newVal);
  virtual void setAtomMass(int index, const float_type& newVal);
  virtual void setAtomFric(int index, const float_type& newVal);
  virtual void setAtomCharge(int index, const float_type& newVal);
  virtual void setSegmentName(int index, const char* newVal);
  virtual void setResidueName(int index, const char* newVal);
  virtual void setResidueNum(int index, const int newVal);
  virtual void setAtomName(int index, const char* newVal);
  virtual void setChemType(int index, const char* newVal);
    
  virtual inline const CDSVec3&        atomPos(int index) const;
  virtual inline const CDSVec3&        atomVel(int index) const;
  virtual inline const float_type& atomMass(int index) const;
  virtual inline const float_type& atomFric(int index) const;
  virtual inline const float_type& atomCharge(int index) const;
  virtual inline const char* segmentName(int index) const;
  virtual inline const char* residueName(int index) const;
  virtual inline int         residueNum(int index) const;
  virtual inline const char* atomName(int index) const;
  virtual inline const char* chemType(int index) const;

  virtual CDSList< int > select(const char* s) const 
  { return subSim()->select(s); }


  float_type calcKE();
  float_type getKE() { return memberData_->ke; }

  struct MemberData {
    pid_t       pid;
    int         barrierCnt;
    int         barrierID;
    float_type  ke;

    // memory allocation comes from shared memory.
    void* operator new(size_t);
    void  operator delete(void*);
    void* operator new(size_t,void* location) { return location; }
    void  operator delete(void*,void*) { }
  };

private:
  MemberData* memberData_;

protected:
  EnsembleSimulation*       ensembleSim_;
  int                       memberIndex_;
  CDS::auto_ptr<CDS::Semaphore>  sem;
  CDSList<int>              ensembleIndices_;
  VectorVec3 atomPosVec_;
  VectorVec3 atomVelVec_;

  friend class EnsembleSimulation;

};


//
// inline methods for EnsembleSimulation
//

inline const Pair& 
EnsembleSimulation::bondPairByID(int bondIndex) const
{ return member()->bondPairByID(bondIndex); }

inline CDSString
EnsembleSimulation::atomString(int index) const
{ return member()->atomString(index); }

inline CDSList< int > 
EnsembleSimulation::select(const char* str) const
{ return member()->select(str); }

inline CDSVector<CDSVec3>
EnsembleSimulation::atomPosArr() const
{ return member()->atomPosArr(); }
inline CDSVector<CDSVec3>
EnsembleSimulation::atomVelArr() const
{ return member()->atomVelArr(); }



inline void
EnsembleSimulation::setAtomPos(int index, const CDSVec3& newVal)
{ member()->setAtomPos(index,newVal); }
inline void
EnsembleSimulation::setAtomVel(int index, const CDSVec3& newVal)
{ member()->setAtomVel(index,newVal); }
inline void
EnsembleSimulation::setAtomMass(int index, const float_type& newVal)
{ member()->setAtomMass(index,newVal); }
inline void 
EnsembleSimulation::setAtomFric(int index, const float_type& newVal)
{ member()->setAtomFric(index,newVal); }
inline void
EnsembleSimulation::setAtomCharge(int index, const float_type& newVal)
{ member()->setAtomCharge(index,newVal); }
inline void
EnsembleSimulation::setSegmentName(int index, const char* newVal)
{ member()->setSegmentName(index,newVal); }
inline void
EnsembleSimulation::setResidueName(int index, const char* newVal)
{ member()->setResidueName(index,newVal); }
inline void
EnsembleSimulation::setResidueNum(int index, const int newVal)
{ member()->setResidueNum(index,newVal); }
inline void
EnsembleSimulation::setAtomName(int index, const char* newVal)
{ member()->setAtomName(index,newVal); }
inline void
EnsembleSimulation::setChemType(int index, const char* newVal)
{ member()->setChemType(index,newVal); }
    
inline const CDSVec3&
EnsembleSimulation::atomPos(int index) const
{ return member()->atomPos(index); }
inline const CDSVec3&
EnsembleSimulation::atomVel(int index) const
{ return member()->atomVel(index); }
inline const float_type&
EnsembleSimulation::atomMass(int index) const
{ return member()->atomMass(index); }
inline const float_type&
EnsembleSimulation::atomFric(int index) const
{ return member()->atomFric(index); }
inline const float_type&
EnsembleSimulation::atomCharge(int index) const
{ return member()->atomCharge(index); }
inline const char*
EnsembleSimulation::segmentName(int index) const
{ return member()->segmentName(index); }
inline const char*
EnsembleSimulation::residueName(int index) const
{ return member()->residueName(index); }
inline int
EnsembleSimulation::residueNum(int index) const
{ return member()->residueNum(index); }
inline const char*
EnsembleSimulation::atomName(int index) const
{ return member()->atomName(index); }
inline const char*
EnsembleSimulation::chemType(int index) const
{ return member()->chemType(index); }

//
// inline methods for EnsembleMemberSimulation
//

inline const CDSVec3 &
EnsembleMemberSimulation::atomPos(int i) const
{
 if ( ensembleSim()->size()==1 )
   return subSim()->atomPos(i);
 
 return atomPosVec_(i);
}

inline const CDSVec3 &
EnsembleMemberSimulation::atomVel(int i) const
{
 if ( ensembleSim()->size()==1 )
   return subSim()->atomVel(i);
 
 return atomVelVec_(i);
}

inline const float_type &
EnsembleMemberSimulation::atomMass(int i) const
{
 return subSim()->atomMass( i );
}

inline const float_type &
EnsembleMemberSimulation::atomFric(int i) const
{
 return subSim()->atomFric( i );
}

inline const char*  
EnsembleMemberSimulation::segmentName(int i) const
{
 return subSim()->segmentName( i );
}

inline const char*
EnsembleMemberSimulation::residueName(int i) const
{
 return subSim()->residueName( i );
}

inline const char*
EnsembleMemberSimulation::chemType(int i) const
{
 return subSim()->chemType( i );
}

inline int
EnsembleMemberSimulation::residueNum(int i) const
{
 return subSim()->residueNum( i );
}

inline const char*
EnsembleMemberSimulation::atomName(int i) const
{
 return subSim()->atomName( i );
}

inline const float_type&
EnsembleMemberSimulation::atomCharge(int i) const
{
 return subSim()->atomCharge( i );
}

//
// pairwise average rmsd between ensemble members
//
class EnsembleRMSD : public AtomSelAction::Base {
  float_type rmsd_;
  CDSMap<int, float_type> byResidue_;
  CDSMap<int, int> count;
  EnsembleSimulation* esim;
public:
  EnsembleRMSD(Simulation* esim=0);

  virtual void run(Simulation* sim,
		   int         atomIndex);

  //overall rmsd
  float_type rmsd() const;  

  //rmsd for each residue
  CDSMap<int, float_type> byResidue() const;

};

#endif /* __ensembleSimulation_hh__ */

#endif
