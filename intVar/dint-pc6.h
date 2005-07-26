#ifndef __dint_pc6__
#define __dint_pc6__

#include "dint-step.h"

class TimeStepAdj;
class VelScale;

/**
 * 6th order predictor-corrector integration.
 */
class PC6 : public Integrator {
public:
    PC6(IVM* ivm) : Integrator(ivm), timeStepAdj(0), velScale(0) {}
    ~PC6();

    void init(const RVec&     pos,
              const RVec&     vel,
              const RVec&     acc);

    virtual void step(double& timeStep);
    static const char* getType() { return "PC6"; }
    bool minimization() const { return 0; }

    void stepUndo();

private:
    RVec         acc;
    RVec         dq3, dq4, dq5;
    RVec         prevAcc;
    TimeStepAdj* timeStepAdj;
    VelScale*    velScale;
    RVec         pos0, vel0, acc0, dq30, dq40, dq50;

    PC6(const PC6&);          //inaccessible
    void operator=(const PC6&);
}; 

#endif /* __dint_pc6__ */

