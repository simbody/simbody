#ifndef __dint_powell__
#define __dint_powell__

#include "dint-step.h"

/**
 * Powell's method minimization.
 */
class Powell : public Integrator {
public:
    Powell(IVM* ivm);
    void init(const RVec& pos,
              const RVec& vel,
              const RVec& acc);
    virtual void step(double& timeStep);
    static const char* getType() { return "Powell"; }
    bool minimization() const { return 1; }

private:
    int     maxcalls;
    double  gradTol;
    double  costTol;
    double  dfpred;
    int     ier;
    int     iterc;
    int     ncalls;
    int     iterfm;
    int     iterrs;
    double  stmin;

    RVec      gcurr;
    double    f;
    RVec      w;
    const int mxfcon;
    const int maxlin;

    double gnew;
    double fmin;
    double gsqrd;
    int    nfopt;
    RVec   wxopt;
    RVec   wgopt;
    double gamden;
    //  double gmin;
    //  double ginit;
    RVec   wrsdg;   // these two
    RVec   wrsdx;   // should be made local
    //
    //  double fx;
    //  double cost;
    //  double costPrev;
    //  RVec   dfx;
    //  RVec   g;
    //  RVec   h;
    //exceptions
    class Exit170 { 
    public: 
        Exit170(int err) : err(err) {} 
        int operator()() { return err; }
    private:
        int err;
    };
    class Exit9000 { 
    public: 
        Exit9000(int err) : err(err) {} 
        int operator()() { return err; } 
    private:
        int err;
    };
    class Exit9005 {};

    double callFG(RVec& pos,
                  RVec& grad );
    void lineMin(const RVec&    w,
                 double&        step,
                 double&        gmin,
                 double&        stepch,
                 double&        sbound,
                 double&        ddspln,
                 const int&     nfbeg,
                 const double&  ginit);
    void getMinVals();
    void saveValsAsOpt(const double& f,
                       const double& sum,
                       const int&    ncalls,
                       const RVec&   pos,
                       const RVec&   gcurr);
    void applyRestart(  RVec&           pos,
                        const double&   beta,
                        const double&   gmin,
                        const double&   ginit,
                        double&         gamden,
                        const RVec&     gcurr,
                        const RVec&     wginit,
                        RVec&           wrsdx,
                        RVec&           wrsdg,
                        RVec&           w,
                        const int&      iterc,
                        int&            iterrs);

    Powell(const Powell&);          //inaccessible
    void operator=(const Powell&);
}; 

#endif /* __dint_powell__ */

