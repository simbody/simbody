#ifndef __dint_step_hh__
#define __dint_step_hh__

/**@file
 * Integration ``step'' routines- for internal dynamics.
 */

#include "dinternal.h"
#include "RVec.h"

class AtomTree;
template<class CHAR> class CDSString;
typedef CDSString<char> String;

/**
 * Abstract base class repesenting "solvers" which are either integrators
 * or minimizers. Concrete numerical methods derive from this class.
 */
class Solver {
public: 
    Solver(IVM* ivm) : ivm(ivm) {}
    virtual ~Solver() {}

    // initialize new solver of given type
    static Solver* create(const String& type,
                          IVM*          ivm);

    virtual void init(const RVec& pos,
                      const RVec& vel,
                      const RVec& acc) = 0;
    virtual bool minimization() const = 0;

    //pos and vel should be taken from member vars at beginning of routine
    //and updated before returning.
    virtual void step(double& timeStep) = 0;
    //  virtual const char* name() = 0;

    virtual RVec& getPos() { return pos; }
    virtual RVec& getVel() { return vel; }
    virtual void setPos(const RVec& p) { pos = p; }
    virtual void setVel(const RVec& v) { vel = v; }

    AtomTree* tree() { return ivm->tree(); }

    //exceptions
    class Finished {
    public: 
        bool ok;
        Finished(bool ok) : ok(ok) {}
    };

protected:
    RVec pos;
    RVec vel;
    IVM*      ivm;
    //  const char* type;
};

#endif /* __dint_step_hh__ */
