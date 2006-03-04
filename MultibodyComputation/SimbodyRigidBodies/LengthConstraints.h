#ifndef LENGTH_CONSTRAINTS_H_
#define LENGTH_CONSTRAINTS_H_

#include "simbody/Simbody.h"
using namespace simtk;

#include "RigidBodyTree.h"

#include <vector>

class RBDistanceConstraint;
class RBDistanceConstraintRuntime;
class RigidBodyTree;

class LengthConstraintsPrivates;
class LengthSet;

class LengthConstraints {
public:
    LengthConstraints(RigidBodyTree&, const double& ctol, int verbose);
    ~LengthConstraints();
    void construct(std::vector<RBDistanceConstraint>&,
                   std::vector<RBDistanceConstraintRuntime>&);

    // Returns true if any change was made in the state.
    bool enforceConfigurationConstraints(SBStateRep&) const;
    bool enforceMotionConstraints(SBStateRep&) const;

    bool calcConstraintForces(const SBStateRep&) const;
    void addInCorrectionForces(const SBStateRep&, SpatialVecList& spatialForces) const;

    void fixVel0(SBStateRep&, Vector&);
    void fixGradient(const SBStateRep&, Vector&);

private:
    //  double tol;
    double bandCut;  //cutoff for calculation of constraint matrix
    int    maxIters;
    int    maxMin;

    RigidBodyTree&             rbTree;
    const int                  verbose;
    LengthConstraintsPrivates* priv;

    friend class LengthSet;
};


#endif /* LENGTH_CONSTRAINTS_H_ */
