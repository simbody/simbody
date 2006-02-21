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

    void enforce(Vector& pos,
                 Vector& vel);

    bool calcConstraintForces(const SBState&) const;
    void addInCorrectionForces(const SBState&, SpatialVecList& spatialForces) const;

    void fixVel0(const SBState&, Vector&);
    void fixGradient(const SBState&, Vector&);

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
