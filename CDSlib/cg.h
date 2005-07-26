
#ifndef __cg_hh__
#define __cg_hh__

//
// conjugate gradient minimizer ala Numerical Recipes Powell algorithm
//
// requires:
//   S    costf(const Vec& x)
//   void gradf(const RVec& x , RVec& grad)
//   bool finished(      int     iter,
//  	           const double& cost,
//		   const RVec&   x,
//		   const RVec&   grad);
//
// CDS - 3/29/00
//

template<class Costf,class Gradf,class Vec,class S,class Condition>
void 
cg(Costf     costf,
   Gradf     gradf,
   Vec&      x,
   S&        cost,
   Condition finished);


template<class Costf,class Gradf,class Vec,class S,class Condition>
void 
cg(Costf     costf,
   Gradf     gradf,
   Vec&      x,
   S&        cost,
   Condition finished)
{
 S fx = costf(x);
 Vec dfx = x; //to properly size gradient.
 gradf(x,dfx);
 dfx *= -1;
 Vec g = dfx;
 Vec h = dfx;
 double stepsize=1.0;
 
 cost = fx;
 for (int iter=1 ; ; iter++) {
   linemin(costf,x,cost,dfx,stepsize);

   if ( finished(iter,cost,x,dfx) )
     return;
   
   fx = costf(x);
   gradf(x,dfx);
   double gg  = abs2(g);
   double dgg = dot( dfx+g , dfx );
   if ( gg == 0.0 ) return;  // unlikely...

   double gam = dgg / gg;
   g = -dfx;
   h *= gam; h += g;
   dfx = h;
 }
} /* cg */

#if 0
// gradient-less version
template<class Costf,class Vec,class S,class Condition>
void 
cg(Costf     costf,
   Vec&      x,
   S&        cost,
   Condition finished)
//p,xi,n,ftol,iter,fret,func)
//float p[],**xi,ftol,*fret,(*func)();
//int n,*iter;
{
//	  int i,ibig,j;
//	  float t,fptt,fp,del;
//	  float *pt,*ptt,*xit,*vector();
//	  void linmin(),nrerror(),free_vector();

//	  pt=vector(1,n);
//	  ptt=vector(1,n);
//	  xit=vector(1,n);

// CDSVector<Vec,1> xi(x.size());
// for (int i=1 ; i<=x.size() ; i++) {
//   xi(i).resize(x.size());
//   xi(i).set(0);
//   xi(i)(i) = 1;
// }
//
// cost = costf(x);

 Vec xt = x;
 Vec xit;
 Vec xtt;
 S fxtt;
 S del;
 for (int iter=1 ; ; iter++) {
   S fx= cost;
   int ibig=0;
   double del=0.0;
   for (int i=1 ; i<=x.size() ; i++) {
     xit = xi(i);
     fxtt= cost;
     linemin(costf,x,cost,xit);
     
     if ( fabs(fxtt-cost) > del ) {
       del = fabs(fxtt-cost);
       ibig=i;
     }
   }
   if ( finished(iter,cost,x,xit) )
     return;

   xtt = 2*x - xt;
   xit = x - xt;
   xt  = x;
   fxtt= costf(xtt);
   if (fxtt < fx) {
     double t=2.0*(fx-2.0*cost+fxtt)*sq(fx-cost-del)-del*sq(fx-fxtt);
     if (t < 0.0) {
       linemin(costf,x,cost,xit);
       xi(ibig) = xit;
     }
   }
 }
} /* cg */
#endif /* 0 */

#endif /* __cg_hh__ */
