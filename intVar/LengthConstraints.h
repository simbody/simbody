#ifndef __dint_loop_hh__
#define __dint_loop_hh__

#include <cdsVector.h>

class IVMAtom;
class AtomClusterNode;
class AtomLoop;
class IVM;
template<class T> class CDSList;

#include <cdsIostream.h>

typedef CDSList<AtomClusterNode*> NodeList;

class LengthConstraintsPrivates;
class LengthSet;

class LengthConstraints {
public:
    LengthConstraints(IVM*);
    ~LengthConstraints();
    void construct(CDSList<AtomLoop>&);
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


#endif /* __dint_loop_hh__ */
