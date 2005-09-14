#ifndef LENGTH_CONSTRAINTS_H_
#define LENGTH_CONSTRAINTS_H_

#include <cdsVector.h>

class RBDistanceConstraint;
class IVM;
template<class T> class CDSList;

#include <cdsIostream.h>

class LengthConstraintsPrivates;
class LengthSet;

class LengthConstraints {
public:
    LengthConstraints(IVM*);
    ~LengthConstraints();
    void construct(CDSList<RBDistanceConstraint>&);

    void enforce(CDSVector<double,1>& pos,
                 CDSVector<double,1>& vel);
    bool fixAccel();
    void fixVel0(CDSVector<double,1>&);
    void fixGradient(CDSVector<double,1>&);

private:
    //  double tol;
    double bandCut;  //cutoff for calculation of constraint matrix
    int    maxIters;
    int    maxMin;

    IVM*                       ivm;
    LengthConstraintsPrivates* priv;

    friend class LengthSet;
};


#endif /* LENGTH_CONSTRAINTS_H_ */
