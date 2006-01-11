#ifndef IVM_LENGTH_CONSTRAINTS_H_
#define IVM_LENGTH_CONSTRAINTS_H_

#include "cdsVector.h"
#include "fixedVector.h"

class IVMDistanceConstraint;
class IVMDistanceConstraintRuntime;
class IVMRigidBodyTree;
template<class T> class CDSList;
typedef FixedVector<double,6>   CDSVec6;
typedef CDSList<CDSVec6> CDSVecVec6;

#include <cdsIostream.h>

class IVMLengthConstraintsPrivates;
class IVMLengthSet;

class IVMLengthConstraints {
public:
    IVMLengthConstraints(IVMRigidBodyTree&, const double& ctol, int verbose);
    ~IVMLengthConstraints();
    void construct(CDSList<IVMDistanceConstraint>&,
                   CDSList<IVMDistanceConstraintRuntime>&);

    void enforce(CDSVector<double,1>& pos,
                 CDSVector<double,1>& vel);

    bool calcConstraintForces() const;
    void addInCorrectionForces(CDSVecVec6& spatialForces) const;

    void fixVel0(CDSVector<double,1>&);
    void fixGradient(CDSVector<double,1>&);

private:
    //  double tol;
    double bandCut;  //cutoff for calculation of constraint matrix
    int    maxIters;
    int    maxMin;

    IVMRigidBodyTree&             rbTree;
    const int                     verbose;
    IVMLengthConstraintsPrivates* priv;

    friend class IVMLengthSet;
};


#endif // IVM_LENGTH_CONSTRAINTS_H_
