#ifndef __dinternal_hh__
#define __dinternal_hh__

//ugly that I have to include these...
#include <cdsList.h>
#include <cdsString.h>
#include <cdsVector.h>
#include <cdsPair.h>

#include "dint-loop.h"
#include "internalDynamics.h"

class IVMAtom;
class AtomTree;
class AT_Build;
class Vec3;
class HingeNode;
class Solver;
class LengthConstraints;

typedef CDSList< IVMAtom* >     AtomList;
typedef CDSVector<double,1>     RVec;   // first element has index 1

typedef int    (*RVecSizeType)(const RVec&);
typedef double (*RVecProdType)(const RVec&, const RVec&);

typedef int    (*VecVec3SizeType)(const CDSVector<Vec3>&);
typedef double (*VecVec3ProdType)(const CDSVector<Vec3>&, 
                                  const CDSVector<Vec3>&);
/**
 * This class owns the atoms and some instructions, builds the
 * multibody system from them and performs all the analyses.
 */
class IVM {
public:
    IVM();
    virtual ~IVM();

    /// Mandatory energy & gradient (forces) routine.
    virtual void calcEnergy() { cout << "you must override this function."; }
    
    /// Kinetic energy and temperature handling (override optional).
    virtual void calcTemperature();

    /// The dimension of an Rvec.
    int rvecSize(const RVec& v) { return rvecSize_(v); }

    /// The inner product of two RVec's.
    double rvecProd(const RVec& v1, const RVec& v2) { return rvecProd_(v1,v2); }

    /// Length squared of an RVec.
    double rvecAbs2(const RVec& v) { return rvecProd(v,v); }

    /// The dimension of a CDSVector<Vec3>.
    int vecVec3Prod(const CDSVector<Vec3>& v) { return vecVec3Size_(v); }

    /// The inner product of two CDSVector<Vec3>'s.
    double vecVec3Prod(const CDSVector<Vec3>& v1,
                       const CDSVector<Vec3>& v2) { return vecVec3Prod_(v1,v2); }

    /// Length squared of a CDSVector<Vec3>.
    double vecVec3Abs2(const CDSVector<Vec3>& v) { return vecVec3Prod(v,v); }

    void initAtoms(const int     natom,
                   const double* massA,
                   const double& kBoltzman);
    void initTopology();
    void initDynamics(bool reuseTopology);

    void groupTorsion();
    void resetCM();
    void printCM();
    void calcY();

    void updateAccel();
    void step(double& timeStep);
    bool minimization() const;

    // get accessors
    const AtomList& getAtoms()  const { return atoms; }
    const AtomTree* tree()      const { return tree_; }
    AtomTree*       tree()            { return tree_; }
    const Solver*   getSolver() const { return solver_; }
    Solver*         getSolver()       { return solver_; }

    int    dof()                  const { return dof_; }
    int    dim()                  const { return dim_; }
    int    verbose()              const { return verbose_; }
    double bathTemp()             const { return bathTemp_; }
    double currentTemp()          const { return currentTemp_; }
    double Etotal()               const { return Etotal_; }
    double Epotential()           const { return Epotential_; }
    double Ekinetic()             const { return Ekinetic_; }
    double Etolerance()           const { return Etolerance_; }
    double Gtolerance()           const { return Gtolerance_; }
    double Ctolerance()           const { return Ctolerance_; }
    double dEpred()               const { return dEpred_; }
    double responseTime()         const { return responseTime_; }
    double frictionCoeff()        const { return frictionCoeff_; }
    double kBoltzmann()           const { return kBoltzmann_; }
    int    maxCalls()             const { return maxCalls_; }
    bool   useLengthConstraints() const { return useLengthConstraints_; }
    bool   adjustTS()             const { return adjustTS_; }
    bool   scaleVel()             const { return scaleVel_; }
    double maxTSFactor()          const { return maxTSFactor_; }
    double maxDeltaE()            const { return maxDeltaE_; }
    double minStepSize()          const { return minStepSize_; }

    // set accessors
    void setVerbose(int v)          { verbose_ = v; }
    void setEtotal(const double& e) { Etotal_ = e; }
    //  void setMinimization(bool m) { minimization_=m; }

    //get descriptive string given atom id
    virtual String idAtom(int id) const;

protected:
    AtomTree*          tree_;
    Solver*            solver_;
    LengthConstraints* lConstraints;

    int dof_;   //number of degrees of freedom
    int dim_;   //number of degrees of freedom+constraints

    int verbose_;
    bool useLengthConstraints_;
    bool adjustTS_;
    bool scaleVel_;
    bool resetCMflag;
    double Etotal_;     //set by calcEnergy
    double Epotential_;
    double Ekinetic_;
    double Etolerance_;
    double Gtolerance_;
    double Ctolerance_;
    double minStepSize_;
    double maxTSFactor_;
    double maxDeltaE_;
    int    maxCalls_;    //maximum allowed calls to energy()
    double kBoltzmann_;
    double bathTemp_;
    double frictionCoeff_;
    double responseTime_;
    double currentTemp_;
    double dEpred_;      //predicted reduction in E - for powell's method

    CDSList<Loop> loops;
    CDSList<Pair> constraintList;

    void groupTorsion(const HingeNode*);
    void initTree();

    AtomList                             atoms;
    CDSList< CDSList<int> >              groupList;
    CDSList<InternalDynamics::HingeSpec> hingeList;
    String                               solverType;
    CDSList<int>                         oldBaseAtoms;

    RVecSizeType rvecSize_;
    RVecProdType rvecProd_;
    VecVec3SizeType vecVec3Size_;
    VecVec3ProdType vecVec3Prod_;
    static int    defaultRVecSize(const RVec&);
    static double defaultRVecProd(const RVec&, const RVec&);
    static int    defaultVecVec3Size(const CDSVector<Vec3>&);
    static double defaultVecVec3Prod(const CDSVector<Vec3>&, 
                                     const CDSVector<Vec3>&);

    friend class AtomTree;
    friend class AT_Build;
    friend class HingeNode;
    template<int DOF> friend class HingeNodeSpec;
};

#endif /* __dinternal_hh__ */
