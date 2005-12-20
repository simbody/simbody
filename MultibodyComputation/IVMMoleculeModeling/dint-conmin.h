#ifndef __dint_conmin__
#define __dint_conmin__

#include "dint-step.h"

/**
 * Conjugate gradient minimization.
 */
class ConMin : public Solver {
public:
    ConMin(IVM*      ivm);

    void init(const RVec& pos,
              const RVec& vel,
              const RVec& acc);
    virtual void step(double& timeStep);
    static const char* getType() { return "ConMin"; }
    bool minimization() const { return 1; }

    double costf(RVec& x);
    void   gradf(RVec& x, RVec& g);
    double costf1D(double& x);

    void linemin(RVec&   p,
                 double& fret,
                 RVec&   xi  ,
                 double& stepsize);
    void mnbrak(double&   ax,
                double&   bx,
                double&   cx,
                double&   fa,
                double&   fb,
                double&   fc);
    double brent(const double& ax,
                 const double& bx,
                 const double& cx,
                 double&       xmin,
                 const double& tol,
                 const int     ITMAX,
                 const double& ZEPS);
private:
    int     maxcalls;  //no yet used
    double  costTol;
    double  sqGradTol;

    //  RVec   x;
    double fx;
    double cost;
    double costPrev;
    RVec   dfx;
    RVec   g;
    RVec   h;
    RVec   pos1D;    // parameters for 1D search
    RVec   dir1D;

    ConMin(const ConMin&);          //inaccessible
    void operator=(const ConMin&);
}; 

#endif /* __dint_conmin__ */

