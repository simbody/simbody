#ifndef __dint_simplex__
#define __dint_simplex__

#include "dint-step.h"

/**
 * Simplex  minimization.
 */
class Simplex : public CDSSolver {
public:
    Simplex(IVM*      ivm);

    void init(const RVec& pos,
              const RVec& vel,
              const RVec& acc);
    virtual void step(double& timeStep);

    static const char* getType() { return "Simplex"; }
    bool minimization() const { return 1; }

    void initSimplex(const double& stepsize);
    double amotry(const int&    ihi,
                  const double& fac);
    double callF(CDSVectorBase<double>& x);
    //  const char* name() { return "Simplex"; }

private:
    typedef CDSVector<double,0> RVec0;

    int                mpts;
    CDSVector<RVec0,0> p;  //points of simplex
    RVec0              y;  //values at points p
    RVec0              psum;
    double             ftol;
    double             lambda;
    int                maxCalls;
    int                ncalls;

    Simplex(const Simplex&);          //inaccessible
    void operator=(const Simplex&);
}; 

#endif /* __dint_simplex__ */

