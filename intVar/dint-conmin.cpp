#include "dint-conmin.h"
#include "dint-step.h"
#include "dinternal.h"
#include <cdsMath.h>

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using namespace CDSMath;
#endif /* USE_CDS_NAMESPACE */


//
// from numerical recipes routine powell
//

using InternalDynamics::printStepDebug;
using InternalDynamics::printStepInfo;

ConMin::ConMin(IVM* ivm)
    : Solver(ivm)
{
}

void
ConMin::init(const RVec& pos,
             const RVec& vel,
             const RVec& acc) 
{
    this->pos = pos;
    maxcalls = ivm->maxCalls();
    costTol = ivm->Etolerance();
    sqGradTol = sq(ivm->Gtolerance());

    fx = costf(this->pos);

    dfx = pos; //to properly size gradient.
    gradf(this->pos,dfx);
    dfx *= -1;
    //   dfx = tree->getInternalForce();
    g = h = dfx;

    cost = fx;
    costPrev = fx;
}

double 
ConMin::costf(RVec& pos) 
{
    RVec vel(pos.size(),0.0);
    tree()->enforceConstraints(pos,vel);;
    tree()->setPosVel(pos,vel);
    ivm->calcEnergy();  //calc energies, derivatives
    if (ivm->verbose()&printStepDebug)
        cout << "costf: " << ivm->Epotential() << '\n';
    return ivm->Epotential();
}

// Presumes that calcEnergy has been called previously with current 
// value of ipos.
void
ConMin::gradf(RVec& x, 
              RVec& g)
{
    g = tree()->getInternalForce();
    // g *= 2.0;
    //   grad *= -1;
}

double 
ConMin::costf1D(double& x)
{
    RVec npos = pos1D + x*dir1D;
    return costf(npos);
}

//
//  Conjugate Gradient Minimizer
//
void
ConMin::step(double& stepSize)
{
    // testGrad(pos);

    cost = fx;
    // x = pos;
    linemin(pos,cost,dfx,stepSize);

    fx = costf(pos);
    if (ivm->verbose()&printStepDebug)
        cout << "step: fx: " << fx << '\n';

    gradf(pos,dfx);
    if (ivm->verbose()&printStepDebug)
        cout << "gg " << abs2(dfx) << '\n';

    if (abs2(dfx) <= sqGradTol ) {
        if (ivm->verbose()&printStepInfo)
            cout << "ConMin::step: exit: gradient tolerance achieved.\n";
        throw Finished(1);
    }
    if ( fabs(cost-costPrev) <= costTol  ) {
        if (ivm->verbose()&printStepInfo)
            cout << "ConMin::step: exit: improvement in cost <= costTol.\n";
        throw Finished(1);
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

void
ConMin::linemin(RVec&   p,
                double& fret,
                RVec&   xi,
                double& stepsize)
{
    double  ax = 0;
    double  bx;
    double  xx = stepsize;
    double  fa, fb, fx;
    pos1D = p;
    dir1D = xi;
    mnbrak(ax,xx,bx,fa,fx,fb);

    // if (ax < 0.0) cout << "linemin: ax<0: " << ax << '\n';
    // if (xx < 0.0) cout << "linemin: xx<0: " << xx << '\n';
    // if (bx < 0.0) cout << "linemin: bx<0: " << bx << '\n';
    // if ( !(ax<xx&&xx<bx) )
    //   cout << "linemin: abscissas out of order: " 
    //	  << ax << ' ' << xx << ' ' << bx << '\n';
    // if ( !( fx<fa && fx<fb) )
    //   cout << "linemin: minimum not bracketed: "
    //	  << fa << ' ' << fx << ' ' << fb << '\n';

    fret = brent(ax,xx,bx,stepsize,costTol,50,1e-8);
    if (ivm->verbose()&printStepDebug)
        cout << "linemin: fret: " << fret << '\n';
    xi *= stepsize;
    p += xi;
}


static const double GOLD   = 1.618034;
static const double GLIMIT = 100.0;
static const double TINY   = 1.0e-20;

inline static
double
SIGN(const double& a,
     const double& b) 
{ return (b > 0.0 ? fabs(a) : -fabs(a)); }

template<class S>
inline static
void
SHFT(   S& a,
        S& b,
        S& c,
        const S& d)
{ a=b; b=c; c=d; }

//template<class S>
//inline static
//void
//swap(      S& a,
//	     S& b)
//{ S tmp=a; a=b; b=tmp; }

// Given function func, and given distinct initial points ax and bx, this
// routine searches in the downhill direction *defined by the function as
// evaluated at the initial points) and returns new points ax, bx, cx which
// bracket a minimum of the function. Also returned are the function values 
// at the three points fa, fb and fc.
void 
ConMin::mnbrak( double& ax,
                double& bx,
                double& cx,
                double& fa,
                double& fb,
                double& fc)
{
    fa=costf1D(ax); //FIX: extra evaluation?
    fb=costf1D(bx);
    if (fb > fa) {
        CDS_NAMESPACE(swap)(ax,bx);
        CDS_NAMESPACE(swap)(fa,fb);
    }
    cx = bx + GOLD*(bx-ax);
    if (ivm->verbose()&printStepDebug)
        cout << "nbrack: initial step: " << cx << ' ';
    fc = costf1D(cx);
    while (fb > fc) {
        double r = (bx-ax) * (fb-fc);
        double q = (bx-cx) * (fb-fa);
        double u = bx - ((bx-cx)*q - (bx-ax)*r) /
                    (2.0*SIGN(max(fabs(q-r),TINY),q-r));
        double ulim = bx + GLIMIT*(cx-bx);
        double fu;
        if ( (bx-u)*(u-cx) > 0.0 ) {
            if (ivm->verbose()&printStepDebug)
                cout << "nbrack: parabolic step: " << u << ' ';
            fu = costf1D(u);
            if (fu < fc) {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            } else if (fu > fb) {
                cx=u;
                fc=fu;
                return;
            }
            u = cx + GOLD*(cx-bx);
            if (ivm->verbose()&printStepDebug)
                cout << "nbrack: min eval: " << u << ' ';
            fu = costf1D(u);
        } else if ( (cx-u)*(u-ulim) > 0.0 ) {
            if (ivm->verbose()&printStepDebug)
                cout << "nbrack: hunt1 eval: " << u << ' ';
            double fu=costf1D(u);
            if (fu < fc) {
                SHFT(bx,cx,u,cx+GOLD*(cx-bx));
                if (ivm->verbose()&printStepDebug)
                    cout << "nbrack: hunt2 eval: " << u << ' ';
                SHFT(fb,fc,fu,costf1D(u));
            }
        } else if ( (u-ulim)*(ulim-cx) >= 0.0 ) {
            u = ulim;
            if (ivm->verbose()&printStepDebug)
                cout << "nbrack: hunt3 eval: " << u << ' ';
            fu = costf1D(u);
        } else {
            u  = cx + GOLD*(cx-bx);
            if (ivm->verbose()&printStepDebug)
                cout << "nbrack: hunt4 eval: " << u << ' ';
            fu = costf1D(u);
        }
        SHFT(ax,bx,cx,u);
        SHFT(fa,fb,fc,fu);
    }
}


static const double CGOLD =  0.3819660;

double
ConMin::brent(  const double& ax,
                const double& bx,
                const double& cx,
                double&       xmin,
                const double& tol,
                const int     ITMAX,
                const double& ZEPS)
{
    double	a=((ax < cx) ? ax : cx);  //a and b bracket minimum w/
    double	b=((ax > cx) ? ax : cx);  //  a<b
    double	x = bx;
    double	w = bx;
    double	v = bx;
    double	fw = costf1D(x);
    double	fv = fw;
    double	fx = fw;

    double	e=0;
    double	d=0;
    for (int iter=1 ; iter<=ITMAX ; iter++) {
        double u,fu;
        double xm=0.5*(a+b);                  //midpoint of a and b
        double tol1 = tol*fabs(x)+ZEPS;
        double tol2= 2 * tol1;
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) { //stop if points close together
            xmin=x;
            return fx;
        }
        if (fabs(e) > tol1) {
            double r=(x-w)*(fx-fv);
            double q=(x-v)*(fx-fw);
            double p=(x-v)*q-(x-w)*r;
            q=2*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            double etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) ) {
                if (ivm->verbose()&printStepDebug)
                    cout << "brent: taking parabolic step: ";
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            } else {
                if (ivm->verbose()&printStepDebug)
                    cout << "brent: taking section step: ";
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        } else {
            if (ivm->verbose()&printStepDebug)
                cout << "brent: taking section step: ";
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=costf1D(u);
        if (fu <= fx) {
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u);
            SHFT(fv,fw,fx,fu);
        } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }
    cerr << "Too many iterations in BRENT\n";
    xmin=x;
    return fx;
}
