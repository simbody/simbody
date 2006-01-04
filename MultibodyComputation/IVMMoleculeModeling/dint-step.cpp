
#include "dint-step.h"
#include "dinternal.h"
#include "AtomTree.h"
#include "dint-powell.h"
#include "dint-conmin.h"
#include "dint-simplex.h"
#include "dint-pc6.h"

#include <cdsMath.h>
#include <linemin.h>
#include <subVector.h>

#include <cdsIomanip.h>

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using CDSMath::sq;
#endif /* USE_CDS_NAMESPACE */

typedef SubVector<RVec>       RSubVec;
typedef SubVector<const RVec> ConstRSubVec;

//
//  Solvers
//

class RungeKutta : public CDSSolver {
public: 
    RungeKutta(IVM*  ivm) : CDSSolver(ivm) {}

    void init(const RVec& pos, const RVec& vel, const RVec& acc) { 
        this->pos = pos; 
        this->vel = vel; 
        ivm->resetCM(); 
    }

    void step(double& timeStep);
    static const char* getType() { return "RungeKutta"; }
    bool minimization() const { return false; }
    RVec deriv(const RVec&);
    //  RVec& getY() { return y; }
};

//
// return the time derivative of v=(pos,vel)
//
RVec
RungeKutta::deriv(const RVec&     v) {
    int dim = v.size()/2;

    RVec pos( ConstRSubVec(v,1    ,dim) );
    RVec vel( ConstRSubVec(v,dim+1,dim) );

    tree()->enforceConstraints(pos,vel);
    tree()->setPosVel(pos,vel);
    ivm->calcEnergy();

    return blockVec(vel, tree()->getAccel());
}

void
RungeKutta::step(double& timeStep)
{
    RVec y = blockVec(getPos(),getVel());

    RVec k1 = timeStep * deriv(y);
    RVec k2 = timeStep * deriv(y+0.5*k1);
    RVec k3 = timeStep * deriv(y+0.5*k2);
    RVec k4 = timeStep * deriv(y+k3);

    y += 1.0/6.0*k1 + 1.0/3.0*k2 + 1.0/3.0*k3 + 1.0/6.0*k4;

    int dim = getPos().size();
    setPos( RVec(RSubVec(y,1    ,dim)) );
    setVel( RVec(RSubVec(y,dim+1,dim)) );

    tree()->enforceConstraints(getPos(),getVel());
    tree()->setPosVel(getPos(),getVel());
    ivm->resetCM();
}

class Gear : public CDSSolver {
    RVec acc;
    RVec acc1, acc2, acc3;
    RVec prevAcc;
    int  stepNum;
    RungeKutta rk;
public:
    //  Gear(RVec& acc) : acc1(acc), acc2(acc), acc3(acc) {} //lazy init
    Gear(IVM*  ivm) : CDSSolver(ivm), rk(ivm) {}

    void init(const RVec& p, const RVec& v, const RVec& a) {
        pos = p;
        vel = v;
        acc = a;
        stepNum = 0;
        rk.init(p,v,a);
    }

    void step(double& timeStep);
    static const char* getType() { return "Gear"; }
    bool minimization() const { return false; }
};

//
// 6th order predictor corrector - internal coordinate version
//
void
Gear::step(double& timeStep)
{
    stepNum++;
    int dim = getPos().size();
    if ( stepNum<4 ) {
        rk.step(timeStep);

        RVec y_  = blockVec( rk.getPos() , rk.getVel() );
        RVec dy_ = rk.deriv(y_);
        RVec acc( RSubVec(dy_,dim+1,dim) );

        acc3 = acc2; acc2 = acc1; acc1 = acc;
        acc1 *= ( 0.5*sq(timeStep) );
    } else {
        //prediction
        vel *= ( timeStep );          //rescale vel/acc
        acc *= ( 0.5*sq(timeStep) );

        prevAcc = acc;

        pos  += vel + 323./180.*acc - 22.0/15.*acc1 + 53./60.*acc2 - 19./90*acc3;
        vel  +=       55./12.  *acc - 59./12. *acc1 + 37./12.*acc2 - 3./4. *acc3;
        acc  +=       3.       *acc - 6.      *acc1 + 4.     *acc2 -        acc3;

        acc3 = acc2; acc2 = acc1; acc1 = prevAcc;

        prevAcc = acc;

        vel *= ( 1.0/timeStep );   //get unscaled velocity
        tree()->enforceConstraints(pos,vel);
        tree()->setPosVel(pos,vel);
        ivm->resetCM();

        ivm->calcEnergy();              //calc energies, derivatives
        acc = tree()->getAccel();
        vel = tree()->getVel(); //may have been tweaked- z.b. normalized. FIX: remove.
        pos = tree()->getPos(); //may have been tweaked- z.b. normalized

        vel *= ( timeStep );          //rescale vel/acc
        acc *= ( 0.5*sq(timeStep) );

        //correction
        RVec dR = acc - prevAcc;

        pos +=   3.0/16.0  * dR;   //the value 3/20 should be used if dim=3*natom
        //   pos +=   3.0/20.0  * dR;   //the value 3/20 should be used if dim=3*natom
        vel += 251.0/360.0 * dR;

        vel *= ( 1.0/timeStep );      //undo vel/acc scaling
        acc *= ( 2.0/sq(timeStep) );
    } 
}
  
//
// Milne predictor-corrector taken from Hamming, Chapter 23
//
class Milne : public CDSSolver {
    //  RVec y0, 
    RVec y1, y2, y3;
    RVec dy0, dy1, dy2;
    double a1, a2, a0;
    double b_, b0, b1, b2;
    double A0, A1, A2, A3;
    double B0, B1, B2;
    int stepNum;
    RungeKutta rk;

public:
//  Milne(int dim) : 
//    y0(dim), y1(dim), y2(dim), y3(dim), 
//    dy0(dim), dy1(dim), dy2(dim),
//    a1(2.0/3.0), a2(1.0/3.0),           //choose for stability
//    A2(0.0), A3(1.0)                 {  //Milne predictor
//  }
    void initConsts() {
        a1 = 2.0/3.0;       //choose for stability
        a2 = 1.0/3.0;                                 
        A2 = 0.0;           //Milne predictor     
        A3 = 1.0;               
        a0 = 1.0 - a1 -a2;
        b_ = (9.0-a1) / 24.0;
        b0 = (19.+13.*a1+8.*a2)/24.0;
        b1 = (-5.+13.*a1+32.*a2)/24.0;
        b2 = (1.0-a1+8.*a2)/24.0;
        A0 = -8.0 - A2 + 8.*A3;
        A1 = 9.0 - 9.0*A3;
        B0 = (17. +   A2 -  9.*A3) / 3.0;
        B1 = (14. + 4*A2 - 18.*A3) / 3.0;
        B2 = (-1. +   A2 +  9.*A3) / 3.0;
    }

    Milne(IVM*  ivm) : CDSSolver(ivm), rk(ivm) {}

    void init(const RVec& p, const RVec& v, const RVec& a) 
    {
        pos = p;
        vel = v;
        dy0 = blockVec(v,a);
        rk.init(p,v,a);
        initConsts();
        //   y1 = y2 = y3 = y0;   //lazy initialization
        //   dy1 = dy2 = dy0;
        stepNum = 0;
    }

    void step(double& timeStep);
    static const char* getType() { return "Milne"; }
    bool minimization() const { return 0; }
};

void
Milne::step(double& timeStep)
{
    stepNum++;

    RVec y0 = blockVec(pos,vel); //previous values
    RVec y_, dy_;                //latest values
    if ( stepNum<4 ) {
        rk.step(timeStep);
        //   y_  = rk.getY();
        pos = rk.getPos();
        vel = rk.getVel();
        y_  = blockVec(pos,vel);
        dy_ = rk.deriv(y_);
    } else {
        //predict
        y_ = A0*y0 + A1*y1 + A2*y2 + A3*y3 +
        timeStep * (B0*dy0 + B1*dy1 + B2*dy2);

        dy_ = rk.deriv(y_);

        //correct
        y_ = a0*y0 + a1*y1 + a2*y2
                + timeStep * (b_*dy_ + b0*dy0 + b1*dy1 + b2*dy2);
        //     dy_ = derivRK(y_);

        int dim = getPos().size();
        pos = RVec( RSubVec(y_,    1,dim) );
        vel = RVec( RSubVec(y_,dim+1,dim) );
    }

    // shift values
    y3 = y2   ; y2 = y1   ; y1 = y0   ; y0 = y_; 
    dy2 = dy1 ; dy1 = dy0 ; dy0 = dy_ ; 
}

class Verlet : public CDSSolver {
    RVec velOldHalf;
public:
    Verlet(IVM*  ivm) : CDSSolver(ivm) {} 
    void init(const RVec& p, const RVec& v, const RVec& a) {
        pos = p;
        vel = v;
        ivm->resetCM();
        velOldHalf = vel;
    }

    void step(double& timeStep);
    static const char* getType() { return "Verlet"; }
    bool minimization() const { return false; }
};

//
// modified Velocity Verlet -- G\"untert, Mumenthaler, W\"uthrich
// J. Mol. Biol. 273, 283 (1997).
//
void
Verlet::step(double& timeStep)
{
    ivm->resetCM();
    ivm->calcEnergy();              //calc energies, derivatives

    double bathFact = 20.0;
    double velScale = sqrt( 1 + (ivm->bathTemp() - ivm->currentTemp())
                        / (bathFact*ivm->currentTemp()) );
    if ( ivm->bathTemp()>=0.0 ) {
        velOldHalf *= velScale;
        vel        *= velScale;
    }
    tree()->setVel(vel);
    tree()->enforceConstraints(pos,vel);

    const RVec a          = tree()->getAccel();
    const RVec velNewHalf = velOldHalf + timeStep * a;
    vel                   = 1.5* velNewHalf - 0.5* velOldHalf;

    pos += timeStep * velNewHalf;
    velOldHalf = velNewHalf;

    tree()->enforceConstraints(pos,vel);
    tree()->setPosVel(pos,vel);
}

//
// from numerical recipes routine powell
//
class MinimizeCG  : public CDSSolver {
    double  costTol;
    double  sqGradTol;

    //  RVec   x;
    double fx;
    double cost;
    double costPrev;
    RVec   dfx;
    RVec   g;
    RVec   h;
public:
    MinimizeCG(IVM* ivm) : CDSSolver(ivm) {}

    void init(const RVec& p, const RVec& v, const RVec& a) {
        pos = p;

        costTol = ivm->Etolerance();
        sqGradTol = sq(ivm->Gtolerance());

        Costf costf(ivm);
        fx = costf(pos);

        dfx = pos; //to properly size gradient.
        fdgradf(pos,dfx);
        dfx *= -1;
        //   dfx = tree->getInternalForce();
        g = h = dfx;

        cost = fx;
        costPrev = fx;
    }

    void step(double& timeStep);
    static const char* getType() { return "MinimizeCG"; }
    bool minimization() const { return true; }

    class Costf {
        IVM*      ivm;
    public: 
        Costf(IVM* ivm) : ivm(ivm) {}
        double operator()(RVec& pos) {
            RVec vel(pos.size(),0.0);
            ivm->tree()->enforceConstraints(pos,vel);;
            ivm->tree()->setPosVel(pos,vel);
            ivm->calcEnergy();  //calc energies, derivatives
            cout << "costf: " << ivm->Epotential() << '\n';
            return ivm->Epotential();
        }
    };

    class Gradf {
        AtomTree* tree;
    public: 
        Gradf(AtomTree* tree) : tree(tree) {}
        void operator()(const RVec& ipos,
        RVec& grad) {
            //presumes that calcEnergy has been called previously with current 
            // value of ipos
            grad = tree->getInternalForce();
            //   grad *= -1;
        }
    };

    void fdgradf(RVec& x, RVec&   grad) {
        //presumes that calcEnergy has been called previously with current 
        // value of ipos
        //   Gradf gradf(tree);
        //   gradf(x,grad); return;
        double eps = 1e-8;
        Costf costf(ivm);
        double f = costf(x);
        for (int i=1 ; i<=grad.size() ; i++) {
            RVec xp = x;
            xp(i) += eps;
            double fp = costf(xp);
            grad(i) = (fp-f) / eps;
        }
        //   grad *= 0.5; //FIX: remove
        //   grad *= -1;
    }

    void testGrad(RVec& x) {
        double tol = 1e-4;
        //   Costf costf(tree);
        //   double f = costf(x);
        RVec grad(x.size());
        Gradf gradf(tree());
        gradf(x,grad);
        RVec fdgrad(x.size());
        fdgradf(x,fdgrad);

        for (int i=1 ; i<=grad.size() ; i++)
            if (fabs(grad(i)-fdgrad(i)) > fdgrad(i)*tol)
                cout << "testGrad: error in gradient: " << setw(3) << i << ' '
                     << grad(i) << ' ' << fdgrad(i) << '\n';
    }
}; 

//
//  Conjugate Gradient Minimizer
//
void
MinimizeCG::step(double& timeStep)
{
    // cg(Costf(),Gradf(),InternalDynamics::pos,cost,
    //    StopCondition(costTol,sqGradTol,numSteps));


    testGrad(pos);

    cost = fx;
    // x = pos;
    linemin(Costf(ivm),pos,cost,dfx,timeStep);

    Costf costf(ivm);
    fx = costf(pos);
    cout << "step: fx: " << fx << '\n';

    fdgradf(pos,dfx);
    cout << "gg " << abs2(dfx) << '\n';

    if (abs2(dfx) <= sqGradTol ) {
        cout << "MinimizeCG::step: exit: gradient tolerance achieved.\n";
        throw Finished(1);
    }
    if ( fabs(cost-costPrev) <= costTol  ) {
        cout << "MinimizeCG::step: exit: improvement in cost <= costTol.\n";
        throw Finished(0);
    }

    costPrev = cost;

    double gg  = abs2(g);
    double dgg = dot( dfx+g , dfx );
    if ( gg == 0.0 ) return;  // unlikely...

    double gam = dgg / gg;
    g = -dfx;
    h *= gam; h += g;
    dfx = h;
}

// CDSSolver factory.
CDSSolver*
CDSSolver::create(const CDSString& type, IVM* ivm)
{
    CDSSolver* ret;
    if ( type.matches(PC6::getType(),true) )
        ret = new PC6(ivm);
    else if ( type.matches(Gear::getType(),true) )
        ret = new Gear(ivm);
    else if ( type.matches(Milne::getType(),true) )
        ret = new Milne(ivm);
    else if ( type.matches(RungeKutta::getType(),true) )
        ret = new RungeKutta(ivm);
    else if ( type.matches(Verlet::getType(),true) )
        ret = new Verlet(ivm);
    else if ( type.matches(MinimizeCG::getType(),true) )
        ret = new MinimizeCG(ivm);
    else if ( type.matches(Powell::getType(),true) )
        ret = new Powell(ivm);
    else if ( type.matches(ConMin::getType(),true) )
        ret = new ConMin(ivm);
    else if ( type.matches(Simplex::getType(),true) )
        ret = new Simplex(ivm);
    else {
        throw IVMInternalDynamics::IVMException("Bad solver specification");
    }
    return ret;
}

