#ifndef LENGTH_CONSTRAINTS_H_
#define LENGTH_CONSTRAINTS_H_

#include "cdsVector.h"
#include "fixedVector.h"

class RBDistanceConstraint;
class RBDistanceConstraintRuntime;
class RigidBodyTree;
template<class T> class CDSList;
typedef FixedVector<double,6>   Vec6;
typedef CDSList<Vec6> VecVec6;

#include <cdsIostream.h>

class LengthConstraintsPrivates;
class LengthSet;

class LengthConstraints {
public:
    LengthConstraints(RigidBodyTree&, const double& ctol, int verbose);
    ~LengthConstraints();
    void construct(CDSList<RBDistanceConstraint>&,
                   CDSList<RBDistanceConstraintRuntime>&);

    void enforce(CDSVector<double,1>& pos,
                 CDSVector<double,1>& vel);

    bool calcConstraintForces() const;
    void addInCorrectionForces(VecVec6& spatialForces) const;

    void fixVel0(CDSVector<double,1>&);
    void fixGradient(CDSVector<double,1>&);

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
