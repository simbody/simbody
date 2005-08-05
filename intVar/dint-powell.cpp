
#include "dint-powell.h"

#include "dint-step.h"

#include "dinternal.h"
#include "AtomTree.h"

#include <cdsMath.h>
#include <cdsString.h>


//Powell's method minimization
// derived from XPLOR code which was taken from an IMSL routine
//

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using namespace CDSMath;
#endif /* USE_CDS_NAMESPACE */


using InternalDynamics::printStepDebug;
using InternalDynamics::printStepInfo;


Powell::Powell(IVM* ivm) : 
  Solver(ivm), mxfcon(2), maxlin(5) 
{}

void
Powell::init(const RVec&     pos,
         const RVec&     vel,
         const RVec&     acc)
{
 this->pos = pos;

 maxcalls = ivm->maxCalls();
 gradTol = ivm->Gtolerance();
 costTol = ivm->Etolerance();
 dfpred = ivm->dEpred();
 ier = iterc = ncalls = iterfm = iterrs = 0;
 gcurr.resize(pos.size());
 f = callFG(this->pos,gcurr);
 w = -gcurr;
 //
 // W(*) IS USED FOR THE SEARCH DIRECTION OF AN ITERATION. WRSDX(*) AND
 // WRSDG(*) CONTAIN THE INFORMATION THAT IS REQUIRED BY THE CONJUGACY
 // CONDITIONS OF THE RESTART PROCEDURE. WGINIT(*) CONTAINS THE GRADIENT
 // AT THE START OF AN ITERATION. WXOPT(*)==(XREF,YREF,ZREF) CONTAINS THE
 // PARAMETERS THAT GIVE THE LEAST CALCULATED VALUE OF F. WGOPT(*)
 // CONTAINS THE GRADIENT VECTOR WHERE F IS LEAST.
 //
 // SET SOME PARAMETERS TO BEGIN THE CALCULATION. ITERC AND NCALLS COUNT
 // THE NUMBER OF ITERATIONS AND CALLS OF ENERGY. ITERFM IS THE NUMBER OF
 // THE MOST RECENT ITERATION THAT DECREASES F.
 //
 //
 // CALL ENERGY. LET THE INITIAL SEARCH DIRECTION BE MINUS THE GRADIENT
 // VECTOR. USUALLY THE PARAMETER ITERRS GIVES THE ITERATION NUMBER OF
 // THE MOST RECENT RESTART, BUT IT IS SET TO ZERO WHEN THE STEEPEST
 // DESCENT DIRECTION IS USED.
 //
 gnew = ivm->rvecProd(w,gcurr);
 double sum  = ivm->rvecAbs2(gcurr);
 if (ivm->verbose()&printStepDebug)
   cout << "Powell::Powell: New search point. GNEW="
    << gnew << " Fnew=" << f << '\n';
 saveValsAsOpt(f,sum,ncalls,pos,gcurr); //35
 
 if ( sum <= gradTol &&          //45
      ivm->verbose()&printStepDebug)
   cout << "Powell::Powell: gradTol achieved at start. "
    << "gsqrd: " << gsqrd << '\n';
 //if ( sum <= gradTol )          //45
 //  throw Finished(1);
 //
 // SET DFPR TO THE ESTIMATE OF THE  REDUCTION IN F GIVEN IN THE
 // ARGUMENT LIST, IN ORDER THAT THE  INITIAL CHANGE TO THE PARAMETERS
 // IS OF A SUITABLE SIZE. THE VALUE  OF STMIN IS USUALLY THE
 // STEP-LENGTH OF THE MOST RECENT  LINE SEARCH THAT GIVES THE LEAST
 // CALCULATED VALUE OF F.
 //
 //CC*** modification.  Original two lines commented out.  
 //CC*** modification by Jiang Jiansheng, York, UK.  12/1/90, ATB. 
 //CC**      DFPR = DFPRED
 //CC**      STMIN = DFPRED/GSQRD
 // dfpred = max( dfpred , gsqrd/pos.size() );
 if ( sum > gradTol )          //45
   stmin = dfpred/gsqrd;
 else
   stmin = dfpred;
} /* Powell::init */

double
Powell::callFG(RVec& pos,
        RVec& grad ) 
  //
  // side effect: pos is modified such that constraints are maintained
  //
{
 ncalls++;
 if ( ncalls>= maxcalls)     //50
   throw Exit9000(131);

 RVec vel(pos.size(),0.0);
 tree()->enforceConstraints(pos,vel);
 tree()->setPosVel(pos,vel);
 ivm->calcEnergy();  //calc energies, derivatives
 grad = tree()->getInternalForce();
 if (ivm->verbose()&printStepDebug) {
   double f = ivm->Epotential();
   double gg = ivm->rvecAbs2(grad);
   cout << "Powell::callFG: Etot=" << f
    << " abs2(g): " << gg << endl;
 }
 return ivm->Epotential();
} /* Powell::callFG */
 

void
Powell::step(double& stepsize)   //FIX: stepsize is not used...
{
 iterc++;             //80
 try {
   if ( ivm->rvecAbs2(gcurr) <= gradTol )
     throw Exit9005();
   //
   // STORE THE INITIAL ENERGY VALUE AND GRADIENT, CALCULATE THE INITIAL
   // DIRECTIONAL DERIVATIVE, AND BRANCH IF ITS VALUE IS NOT NEGATIVE. SET
   // SBOUND TO MINUS ONE TO INDICATE THAT A BOUND ON THE STEP IS NOT KNOWN
   // YET, AND SET NFBEG TO THE CURRENT VALUE OF NCALLS. THE PARAMETER
   // IRETRY SHOWS THE NUMBER OF ATTEMPTS AT SATISFYING THE BETA CONDITION.
   //
   double finit  = f;
   RVec   wginit = gcurr;
   double ginit  = ivm->rvecProd(w,gcurr);
   if ( ginit>=0.0 ) 
     throw Exit170(130);
   double gmin = ginit;
   double sbound = -1.0;
   int    nfbeg = ncalls;
   //
   // SET STEPCH SO THAT THE INITIAL STEP-LENGTH IS CONSISTENT WITH THE
   // PREDICTED REDUCTION IN F, SUBJECT TO THE CONDITION THAT IT DOES NOT
   // EXCEED THE STEP-LENGTH OF THE PREVIOUS ITERATION. LET STMIN BE THE
   // STEP TO THE LEAST CALCULATED VALUE OF F.
   //
   double stepch = min(stmin,fabs(dfpred/ginit));
   stmin = 0.0;
   //
   // CALL ENERGY AT THE VALUE OF X THAT IS DEFINED BY THE NEW CHANGE TO
   // THE STEP-LENGTH, AND LET THE NEW STEP-LENGTH BE STEP.
   //
   double stepsize = stmin+stepch;  //90
   pos = wxopt + stepch*w;
   double ddspln=0.;
   int    iretry = -1;
   RVec diff = pos-wxopt;
   if ( ivm->rvecAbs2(diff) > 0.0 ) {
     lineMin(w,stepsize,gmin,stepch,sbound,ddspln,nfbeg,ginit);
   } else 
     if (ncalls>nfbeg+1 ||
     fabs(gmin/ginit)>0.2) {
       if (ivm->verbose()&printStepDebug)
     cout << "Powell::step: abandoned: pos=wxopt.\n";
       throw Exit170(129);
       //       iretry=1;
       //       break;
     }
   double sum = 0.0;
   double beta = 0.0;
   while (1) {
     getMinVals();            //170
     //
     // CALCULATE THE VALUE OF BETA THAT OCCURS IN THE NEW SEARCH
     // DIRECTION.
     //
     sum = ivm->rvecProd(gcurr,wginit);  //125
     beta = (gsqrd-sum)/(gmin-ginit);
     //
     // TEST THAT THE NEW SEARCH DIRECTION CAN BE MADE DOWNHILL. IF IT
     // CANNOT, THEN MAKE ONE ATTEMPT TO IMPROVE THE ACCURACY OF THE LINE
     // SEARCH.
     //
     if ( fabs(beta*gmin)<=0.2*gsqrd )
       break;                   //exit to 135
     iretry++;
     if ( iretry>0 )
       break;                   //exit to 135
     if ( ncalls>=nfopt+maxlin ) {    //something wrong with line search
       iretry=1;                      //use steepest descent dir
       break;
     }
     if ( ncalls>=nfopt+maxlin ) { //110
       if (ivm->verbose()&printStepDebug)
     cout << "Powell::step: beta loop: ncalls>=nfopt+maxlin.\n";
       throw Exit170(129);      //change 1
       //       break;
     }
     stepch = 0.5*(sbound-stmin); //120
     if (sbound < -0.5) 
       stepch = 9.0*stmin;
     double gspln = gmin+stepch*ddspln;
     if (gmin*gspln<0.0)
       stepch *= gmin/(gmin-gspln);
     stepsize = stmin+stepch;         //90
     pos = wxopt + stepch*w;
     if ( ivm->rvecAbs2(pos - wxopt) > 0.0 )
       lineMin(w,stepsize,gmin,stepch,sbound,ddspln,nfbeg,ginit);
     else
       if ( ncalls>nfbeg+1      ||
        fabs(gmin/ginit)>0.2 ) {
     if (ivm->verbose()&printStepDebug)
       cout << "Powell::step: lineMin abandoned.\n";
     Exit170(129);
       }
   }
   //
   // APPLY THE TEST THAT DEPENDS ON THE PARAMETER MXFCON. SET DFPR TO THE
   // PREDICTED REDUCTION IN F ON THE NEXT ITERATION.
   //
   if (f<finit)                  //135
     iterfm = iterc;
   if (iterc >= iterfm+mxfcon)
     throw Exit9000(132);
   dfpred = stmin*ginit;           //140

   //
   // BRANCH IF A RESTART PROCEDURE IS REQUIRED DUE TO THE ITERATION NUMBER
   // OR DUE TO THE SCALAR PRODUCT OF CONSECUTIVE GRADIENTS.
   //
   if ( iretry>0 ) {
     // use stepest descent direction
     w = -gcurr;                 //10
     tree()->enforceConstraints(pos,w);  //FIX: should already be enforced
     iterrs = 0;
   } else if (iterrs==0                ||
          iterc-iterrs>=ivm->rvecSize(pos) ||
          fabs(sum)>=0.2*gsqrd      ) {
     applyRestart(pos,
          beta,gmin,ginit,gamden,
          gcurr,wginit,wrsdx,wrsdg,w,
          iterc,iterrs);
   } else {
     //
     // CALCULATE THE VALUE OF GAMA THAT OCCURS IN THE NEW SEARCH DIRECTION,
     // AND SET SUM TO A SCALAR PRODUCT FOR THE TEST BELOW. THE VALUE OF
     // GAMDEN IS SET BY THE RESTART PROCEDURE.
     //
     double gamma = ivm->rvecProd(gcurr,wrsdg) / gamden;
     double sum   = ivm->rvecProd(gcurr,wrsdx);
     //
     // RESTART IF THE NEW SEARCH DIRECTION IS NOT SUFFICIENTLY DOWNHILL.
     //
     if ( fabs(beta*gmin + gamma*sum) >= 0.2*gsqrd )
       applyRestart(pos,
            beta,gmin,ginit,gamden,
            gcurr,wginit,wrsdx,wrsdg,w,
            iterc,iterrs);
     else {
       //
       // CALCULATE THE NEW SEARCH DIRECTION.
       //
       if (ivm->verbose()&printStepDebug)
     cout << "Powell::step: New search direction. beta="
          << beta << " GAMA=" << gamma << '\n';
       w = -gcurr + beta*w + gamma*wrsdx;
       tree()->enforceConstraints(pos,w);
     }
   }
 }
 catch ( InternalDynamics::Exception e) {
   if ( String(e.mess).contains("LengthSet::enforce") ) {
     getMinVals();
     if (ivm->verbose()&printStepInfo) 
       cout << "Powell::step: step terminated because of unmet constraint.\n";
   } else 
     throw;
 }
  
       
 catch ( Exit170 error ) {
   getMinVals();
   RVec vel(pos.size(),0.0);
   tree()->enforceConstraints(pos,vel);
   tree()->setPosVel(pos,vel);
   if (ivm->verbose()&printStepInfo) {
     cout << "Powell::step: irregular exit: ";
     switch ( error() ) {
       case 129: 
     cout << "Line search abandoned: gradient may be incorrect\n";
     break;
       case 130: 
     cout << "Search direction uphill\n";
     break;
     }
   }
   throw Finished(0);
 }
 catch ( Exit9000 error ) {
   getMinVals();
   RVec vel(pos.size(),0.0);
   tree()->enforceConstraints(pos,vel);
   tree()->setPosVel(pos,vel);
   if (ivm->verbose()&printStepInfo) {
     cout << "Powell::step: irregular exit: ";
     switch ( error() ) {
       case 133: 
     cout << "energy reduction is less than costTol\n";
     break;
       case 132: 
     cout << "Two consecutive iterations failed to reduce E\n";
     break;
       case 131: 
     cout << "max number of energy evaluations reached\n";
     break;
     }
   }
   throw Finished(0);
 }
 catch ( Exit9005 ) {
   if (ivm->verbose()&printStepInfo)
     cout << "Powell::step: normal exit.\n";
   throw Finished(1);
 }

} /* Powell::step */

void
Powell::applyRestart(      RVec&   pos,   
              const double& beta,
              const double& gmin,
              const double& ginit,
                double& gamden,
              const RVec&   gcurr,
              const RVec&   wginit,
                RVec&   wrsdx,
                RVec&   wrsdg,
                RVec&   w,
              const int&    iterc,
                int&    iterrs)
  //
  // label 155 from powell.s
  //
{
 gamden = gmin-ginit;
 if (ivm->verbose()&printStepDebug)
   cout << "Powell::applyRestart: Applying restart procedure. BETA="
    << beta << '\n';
 wrsdx = w;
 wrsdg = gcurr-wginit;
 w     = -gcurr + beta*w;
 tree()->enforceConstraints(pos,w); //should not change pos
 //callFG(); //FIX: necessary?
 iterrs = iterc;
} /* Powell::applyRestart */

void
Powell::getMinVals()
{
 if (ncalls != nfopt) {             // 170
   f = fmin;
   pos = wxopt;
   gcurr = wgopt;
   if (ivm->verbose()&printStepDebug)
     cout << "Powell::saveMinVals: Current point set to optimal point.\n"
      << "\tf: " << f << " abs2(g): " << ivm->rvecAbs2(gcurr) << endl;
 }
} /* Powell::getMinVals */

void
Powell::saveValsAsOpt(const double& f,
               const double& sum,
               const int&    ncalls,
               const RVec&   pos,
               const RVec&   gcurr)
  //
  // corresponds to label 35
  //
{
 fmin = f;
 gsqrd = sum;
 nfopt = ncalls;
 wxopt = pos;
 wgopt = gcurr;
 if (ivm->verbose()&printStepDebug)
   cout << "Powell::saveValsAsOpt: least energy point set to current point.\n";
} /* Powell::saveValsAsOpt */

void
Powell::lineMin(const RVec&   w,
               double& stepsize,
               double& gmin,
               double& stepch,
               double& sbound,
               double& ddspln,
         const int&    nfbeg,
         const double& ginit)
  //
  // find minimum along the w direction
  //
  // ( program block beginning at 5 in XPLOR's powell.s )
  //
{
 int startCall = ncalls;
 while (1) {
   f = callFG(pos,gcurr);
   gnew = ivm->rvecProd(w,gcurr);
   double sum = ivm->rvecAbs2(gcurr);
   if (ivm->verbose()&printStepDebug)
     cout << " Powell::lineMin: New search point. GNEW=" << gnew
      << " Fnew=" << f << " gg=" << sum << '\n';
   double fch = f-fmin;
   //
   // STORE THE VALUES OF XCURR, F AND GCURR, IF  THEY ARE THE BEST THAT
   // HAVE BEEN CALCULATED SO FAR, AND NOTE GCURR  SQUARED AND THE VALUE OF
   // NCALLS. TEST FOR CONVERGENCE.
   //
   if ( fabs(fch) < costTol )   //new
     throw Exit9000(133);
   if (fch < 0.0) {
     saveValsAsOpt(f,sum,ncalls,pos,gcurr);     //35
     if ( sum <= gradTol ) 
       throw Exit9005();       //45
   } else if (fch == 0.0) {
     //     if ( gnew/gmin >= -1.0 )  //original logic
     if ( gnew/gmin < 0.0 )           //gradients of diff. sign: diff. wells.
       saveValsAsOpt(f,sum,ncalls,pos,gcurr);   //35
     if ( sum <= gradTol ) 
       throw Exit9005();       //45
   }
   if ( ncalls>= maxcalls)     //50
     throw Exit9000(131);
   
   //
   // Let spln be the quadratic spline that interpolates the calculated
   // energy values and directional derivatives at the points stmin and
   // step of the line search, where the knot of the spline is at
   // 0.5*(stmin+step). revise stmin, gmin and sbound, and set ddspln to
   // the second derivative of spln at the new stmin. However, if fch is
   // zero, it is assumed that the maximum accuracy is almost achieved, so
   // ddspln is calculated using only the change in the gradient.
   //
   double tmp    = 2*fch/stepch-gnew-gmin;  //100
   ddspln = (gnew-gmin)/stepch;
   
   if ( ncalls > nfopt )  // most recent pos didn't have lowest E
     sbound = stepsize;
   else {                 // ncalls == nfopt: most recent pos at lowest E
     if ( gmin*gnew<=0.0 ) 
       sbound = stmin;
     stmin = stepsize;
     gmin = gnew;
     stepch *= -1;
   }
   if (fch != 0.0)                    //105
     ddspln += 2*tmp/stepch;
   //
   // TEST FOR CONVERGENCE OF THE LINE SEARCH, BUT FORCE AT LEAST TWO STEPS
   // TO BE TAKEN IN ORDER NOT TO LOSE QUADRATIC TERMINATION.
   //
   if ( gmin == 0.0 )
     break;
   if (ncalls > nfbeg+1 ) {
     if ( fabs(gmin/ginit)<= 0.2 )
       break;                       //exit to 170
     if (ncalls>=nfopt+maxlin) {     //110
       if (ivm->verbose()&printStepDebug)
     cout << "Powell::lineMin: ncalls>=nfopt+maxlin. "
          << "Line min abandoned.\n";
       throw Exit170(129);      //change 2
       //break;     //abandon line search
     }
   }
   //
   // SET STEPCH TO THE GREATEST CHANGE TO THE CURRENT VALUE OF STMIN THAT
   // IS ALLOWED BY THE BOUND ON THE LINE SEARCH. SET GSPLN TO THE GRADIENT
   // OF THE QUADRATIC SPLINE AT (STMIN+STEPCH). HENCE CALCULATE THE VALUE
   // OF STEPCH THAT MINIMIZES THE SPLINE FUNCTION, AND THEN OBTAIN THE NEW
   // FUNCTION AND GRADIENT VECTOR, FOR THIS VALUE OF THE CHANGE TO THE
   // STEP-LENGTH.
   //
   stepch = (sbound<-0.5) ? 9.0*stmin : 0.5*(sbound-stmin);   //120
   double gspln = gmin+stepch*ddspln;
   if (gmin*gspln<0.0)
     stepch *= gmin/(gmin-gspln);
   //
   // CALL ENERGY AT THE VALUE OF X THAT IS DEFINED BY THE NEW CHANGE TO
   // THE STEP-LENGTH, AND LET THE NEW STEP-LENGTH BE STEP. 
   //
   stepsize = stmin+stepch;              //90
   pos = wxopt + stepch*w;
   if ( ivm->rvecAbs2( pos-wxopt ) == 0.0 ) {
     if (ivm->verbose()&printStepDebug)
       cout << "Powell::lineMin: pos=wxopt. Line min abandoned.\n";
     if (ncalls>nfbeg+1     ||
     fabs(gmin/ginit)>0.2)
       throw Exit170(129);
     else
       break;
   }
 }
 if ( nfopt <= startCall ) //made no improvement
   throw Exit170(129);
}

