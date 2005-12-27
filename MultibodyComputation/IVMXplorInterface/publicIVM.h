#ifndef __publicIVM_hh__
#define __publicIVM_hh__

// public interface to the IVM class
// appropriate for scripting interfaces (and SWIG understands this file).
// 


#include <sthead.h>
#include <cdsExcept.h>
#include <simulation.h>
#include <enumNameMap.h>
#include "dinternal.h"
#include "internalDynamics.h"
#include "publicNode.h"
#include "potList.h"
#include <modified.h>


using InternalDynamics::HingeSpec;

class PublicNode;


//FIX: tree is not defined until init is called. 
// members such as setPos will not work in this state (segfault). Should
// another python object be created?

class PublicIVM : public IVM, public ModifiedBase {
public:
    PublicIVM(Simulation* sim)
      : IVM(), finalTime(0.), steps(0), eCount_(0), sim(sim)
    {  registerTo(sim); }
    //  ~PublicIVM() { unRegister(sim); }

    void init();
    void initDynamics(bool reuseTopology);

    // take dynamics or minimization step
    //   TRUE is returned if final condition is met
    //   if stepsize is variable, stepsize argument is modified
    bool step(float_type& stepsize);

    // group atoms in a manner appropriate for torsion-angle dynamics
    //   obeys pre-existing groupList, hingeList, bondExclude and constraintList
    void groupTorsion();


    //energy and derivatives calculated and terms placed in energyTerms
    //  also, eCount is incremented
    void calcEnergy();
    void calcTemperature();

    //  CDSList<EnergyReport> energyTerms();
    void eCountReset() { eCount_=0; }

    //zero center of mass linear and angular momentum
    void resetCM() { IVM::resetCMflag=1; }

    //given an atom index, return a descriptive string
    virtual CDSString idAtom(int id) const;

    //bonds to ignore when considering molecular topology
    void setBondExclude(const CDSList<Pair>& list ) { bondExclude_=list; }

    //list of list of atom indices- grouped atoms fixed w/ respect to each other
    void setGroupList(const CDSList< CDSList<int> >& list);

    //bonds to constrain (if useLengthConstraints)
    void setConstraintList(const CDSList<Pair>& list );

    //atom(s) to use as base atoms in the tree topology
    void setOldBaseAtoms(const CDSList<int>& list );

    //values of the internal coordinates (NOT atom positions)
    void setPos(const CDSVector<double>& pos);

    // get list with info on individual nodes
    CDSList<PublicNode> nodeList();

    //specification of hinge type assigned to atoms
    //  note that the specified atoms are not grouped together in fixed
    //  nodes: use groupList for that purpose.
    void setHingeList(const CDSList<HingeSpec>& list );

    //set internal variable velocities to be as close as possible to
    // the current Cartesian counterparts
    void velFromCartesian();

    // associated get-accessors
    CDSList<Pair>            bondExclude() const { return bondExclude_; }
    CDSList< CDSList<int> >  groupList();
    CDSList<Pair>            constraintList();
    CDSList<int>             oldBaseAtoms();
    CDSList<HingeSpec>       hingeList();
    CDSVector<double>        pos() const;

    //simple get-accessors
    int        dof()                  const { return dof_; }
    int        dim()                  const { return dim_; }
    int        verbose()              const { return verbose_;}
    bool       minimization()         const { return IVM::minimization(); }
    int        eCount()               const { return eCount_; }
    float_type bathTemp()             const { return bathTemp_;}
    float_type currentTemp()          const { return currentTemp_; }
    float_type Etotal()               const { return Etotal_; }
    float_type Epotential()           const { return Epotential_; }
    float_type Ekinetic()             const { return Ekinetic_; }
    float_type gradMagnitude()        const { return gradMagnitude_; }
    float_type eTolerance()           const { return Etolerance_; }
    float_type gTolerance()           const { return Gtolerance_; }
    float_type cTolerance()           const { return Ctolerance_; }
    float_type dEpred()               const { return dEpred_; }
    float_type responseTime()         const { return responseTime_; }
    float_type frictionCoeff()        const { return frictionCoeff_; }
    float_type kBoltzmann()           const { return kBoltzmann_; }
    int        maxCalls()             const { return maxCalls_; }
    bool       constrainLengths()     const { return useLengthConstraints_; }
    bool       adjustStepsize()       const { return adjustTS_; }
    bool       scaleVel()             const { return scaleVel_; }
    float_type maxTSFactor()          const { return maxTSFactor_; }
    float_type maxDeltaE()            const { return maxDeltaE_; }
    float_type minStepSize()          const { return IVM::minStepSize(); }
    const char* stepType()            const { return IVM::solverType; }
    const PotList& potList()          const { return potList_; }
    PotList&   potList()                    { return potList_; }

    // simple set-accessors
    void setVerbose(const int v)               { verbose_=v;}
    void setConstrainLengths(bool i)           { useLengthConstraints_ = i; }
    void setStepType(const char* v)            { IVM::solverType=v; }
    void setBathTemp(const float_type& v)      { bathTemp_=v; }
    void setETolerance(const float_type& v)    { Etolerance_=v; }
    void setGTolerance(const float_type& v)    { Gtolerance_=v; }
    void setCTolerance(const float_type& v)    { Ctolerance_=v; }
    void setDEpred(const float_type& v)        { dEpred_=v; }
    void setResponseTime(const float_type& v)  { responseTime_=v; }
    void setFrictionCoeff(const float_type& v) { frictionCoeff_=v; }
    void setMaxCalls(const int i)              { maxCalls_=i; }
    void setAdjustStepsize(const bool b)       { adjustTS_=b; }
    void setScaleVel(const bool b)             { scaleVel_=b; }
    void setMaxTSFactor(const float_type& v)   { maxTSFactor_=v; }
    void setMaxDeltaE(const float_type& v)     { maxDeltaE_=v; }
    void setMinStepSize(const float_type& v)   { minStepSize_=v; }
    void setPotList(const PotList& p)          { potList_ = p; }

    // insert constants from InternalDynamics namespace
    enum VerboseFlags {
        printCoords          = InternalDynamics::printCoords,          
        printResetCM         = InternalDynamics::printResetCM,         
        printVelFromCartCost = InternalDynamics::printVelFromCartCost, 
        printTemperature     = InternalDynamics::printTemperature,     
        printEnergy          = InternalDynamics::printEnergy,          
        printCMVel           = InternalDynamics::printCMVel,           
        printNodeForce       = InternalDynamics::printNodeForce,       
        printNodePos         = InternalDynamics::printNodePos,         
        printNodeTheta       = InternalDynamics::printNodeTheta,       
        printStepDebug       = InternalDynamics::printStepDebug,       
        printStepInfo        = InternalDynamics::printStepInfo,        
        printNodeDef         = InternalDynamics::printNodeDef,         
        printLoopDebug       = InternalDynamics::printLoopDebug,       
        printLoopInfo        = InternalDynamics::printLoopInfo
    };  

    //  ostream& print(ostream& ostr);

    //callback setup
    void setRVecSize(RVecSizeType cb) { rvecSize_ = cb; }
    void setRVecProd(RVecProdType cb) { rvecProd_ = cb; }
    void setVecVec3Size(VecVec3SizeType cb) { vecVec3Size_ = cb; }
    void setVecVec3Prod(VecVec3ProdType cb) { vecVec3Prod_ = cb; }

private:
    float_type      finalTime;
    int             steps;
    int             eCount_;        // number times calcEnergy is called in current step
    float_type      gradMagnitude_; // magnitude of gradient
    CDSList<Pair>   bondExclude_;
    PotList         potList_;
    Simulation*     sim;
    CDSVector<CDSVec3> shadowPos;      //vectors shadowing simulation values
    CDSVector<CDSVec3> shadowVel;

    //  sync members for private interface
    void syncPos();
    void syncVel();
    void syncPosFrom();
    void syncVelFrom();

    void updateValues();
};

namespace EnumNamespace {
    extern EnumNameMap VerboseFlags[];
};

#endif /* __publicIVM_hh__ */
