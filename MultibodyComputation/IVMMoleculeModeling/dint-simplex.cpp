//
// FIX: needs to be tested. Is the expression for lambda reasonable?
//

#include "dint-simplex.h"

#include "dint-step.h"

#include "dinternal.h"
#include "AtomTree.h"

#include <cdsMath.h>

#include <cdsVector.h>


const static double ALPHA = 1.0;
const static double BETA  = 0.5;
const static double GAMMA = 2.0;



using InternalDynamics::printStepDebug;
using InternalDynamics::printStepInfo;


Simplex::Simplex(IVM* ivm) : CDSSolver(ivm) {}

void
Simplex::init(const RVec& pos,
	      const RVec& vel,
	      const RVec& acc)
{
 this->pos = pos;
 mpts = pos.size()+1;
 p.resize(mpts);
 y.resize(mpts);
 psum.resize(pos.size());
 psum.set(0.0);
 ftol = ivm->Etolerance();
 maxCalls = ivm->maxCalls();
 ncalls = 0;

 p(0)=pos;
 y(0) = callF(p(0));
 lambda = 0.;
} /* Simplex::Simplex */

void
Simplex::initSimplex(const double& stepsize)
{
 lambda = stepsize;
// if (pos.size()>1) {
//   lambda = 0;
//   for (int i=1 ; i<pos.size() ; i++) 
//     lambda += sq(pos(i)-pos(i-1));
//   lambda = sqrt(lambda) / (pos.size()-1);
// }
 if (ivm->verbose()&printStepDebug)
   cout << "Simplex::Simplex: pos = " << pos 
	<< "\n\tlambda initialized to " << lambda << '\n';
 for (int i=0 ; i<pos.size() ; i++) {
   RVec0 ei(pos.size(),0.0);
   ei(i) = 1.0;
   p(i+1) = pos + lambda * ei;
   y(i+1) = callF(p(i+1));
 }
 for (int j=0 ; j<psum.size() ; j++) { 
   for (int i=0 ; i<mpts ; i++) 
     psum(j) += p(i)(j); 
 } 
} /* Simplex::init */

void
Simplex::step(double& stepsize) 
{
 if ( lambda == 0. )
   initSimplex(stepsize);

 int ilo = 0;                /* determine which point is the highest (worst),*/
 int inhi=0;
 int ihi = y(0)>y(1) ? (inhi=1,0) : (inhi=0,1);   /*next-highest, and lowest */
 for (int i=0 ; i<mpts ; i++) {            /*best-by looping over the points */
   if (y(i) < y(ilo)) ilo = i;
   if (y(i) > y(ihi)) {
     inhi = ihi;
     ihi = i;
   } else if (y(i) > y(inhi)) if (i != ihi) inhi = i;
 }

 /*print out value at the lowest point*/
 if (ivm->verbose()&printStepDebug)
   cout << "Simplex::step: Value at lowest point: " << y(ilo) << '\n';

 double rtol = 2.0 * fabs(y(ihi)-y(ilo)) / (fabs(y(ihi))+fabs(y(ilo)));
 /*Compute the fractional range from highest to lowest and return if
   satisfactory*/
 if (rtol < ftol) {
   if (ivm->verbose()&printStepInfo)
     cout << "Simplex::step: normal exit.\n";
   pos = p(ilo);
   throw Finished(1);
 }
 if (ncalls >= maxCalls) {
   if (ivm->verbose()&printStepInfo)
     cout << "Simplex::step: abnormal exit: maxCalls exceeded.\n";
   pos = p(ilo);
   throw Finished(0);
 }
 /*Begin a new iteration. First extrapolate by a fcator ALPHA through
   the face of the simplex across from the high point, i.e., reflect the
   simplex from the high point.*/  
 double ytry = amotry(ihi,-ALPHA);
 if (ytry <= y(ilo))
   /*Gives a result better than the best point, so try an additional
     extrapolation by a factor GAMMA.*/
   ytry = amotry(ihi,GAMMA);
 else if (ytry >=y(inhi)) {
   /* The reflected point is worse than the second-highest, so look
      for an intermediate lower point, i.e., do a one-dimensional contraction.*/
   double ysave = y(ihi);
   ytry = amotry(ihi,BETA);
   if (ytry >= ysave) { 
     /*Can't seem to get rid of that high point. Better contract
       around the lowest (best) point.*/
     if (ivm->verbose()&printStepDebug)
       cout << "Simplex::step: contracting.\n";
     for (int i=0;i<mpts;i++) {
       if (i!= ilo) {
	 p(i) = 0.5 * (p(i)+p(ilo));
	 y(i) = callF(psum);
       }
     }
     psum.set( 0.0 );
     for (int j=0 ; j<psum.size() ; j++) 
       for (int i=0 ; i<mpts ; i++) 
	 psum(j) += p(i)(j); 
   }
 }
 pos = p(ilo);
} /* Simplex::step */

double
Simplex::callF(CDSVectorBase<double>& x) 
{
 RVec pos(x);
 ncalls++;

 RVec vel(pos.size(),0.0);
 tree()->enforceConstraints(pos,vel);
 tree()->setPosVel(pos,vel);
 ivm->calcEnergy();  //calc energies, derivatives
 if (ivm->verbose()&printStepDebug)
   cout << "Simplex::callF: Etot=" << ivm->Epotential() << '\n';
 return ivm->Epotential();
} /* Simplex::callF */

double
Simplex::amotry(const int&    ihi,
		 const double& fac)
  //Extrapolate by a factor fac through the face of the simplex across
  //    from the high point, tries it and replaces the high point if the new
  //    point is better
{
 double fac1 = (1.0-fac) / pos.size();
 double fac2 = fac1-fac;
 RVec0  ptry = fac1*psum-fac2*p(ihi);
 double ytry = callF(ptry);    /*Evaluate the function at the trial point*/
 if (ytry < y(ihi)) { 
   /*If its better than the highest, then replace the highest.*/
   y(ihi) = ytry;
   psum += ptry - p(ihi);
   p(ihi) = ptry;
 }
 return ytry;
} /* Simplex::amotry */
