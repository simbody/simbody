#ifndef __dint_loop_hh__
#define __dint_loop_hh__

#include <cdsVector.h>

class IVMAtom;
class HingeNode;
class Loop;
class IVM;
template<class T> class CDSList;

#include <cdsIostream.h>

typedef CDSList<HingeNode*> NodeList;

class LengthConstraintsPrivates;
class LengthSet;

class LengthConstraints {
public:
    LengthConstraints(IVM*);
    ~LengthConstraints();
    void construct(CDSList<Loop>&);
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

class Loop {
public:
    Loop() {}
    Loop(IVMAtom* baseAtom,
         IVMAtom* tipAtom);

protected:
    IVMAtom* tip1;
    IVMAtom* tip2;

    friend class LoopWNodes;
    friend class LengthSet;
    friend void LengthConstraints::construct(CDSList<Loop>& iloops);
    friend ostream& operator<<(ostream& os, const LengthSet& s);
};

#endif /* __dint_loop_hh__ */
