
#ifndef __mnbrak_hh__
#define __mnbrak_hh__

template<class S,class T, class Func>
void 
mnbrak(S&   ax,
       S&   bx,
       S&   cx,
       T&   fa,
       T&   fb,
       T&   fc,
       Func func);
//#include <math.h>

#ifdef MAX
#undef MAX
#endif


template<class T>
inline static
T 
MAX(const T& a,
    const T& b)
{ return (a > b ? a : b); }

template<class S,class T>
inline static
const T
SIGN(const T& a,
     const S& b) 
{ return ((b) > 0.0 ? fabs(a) : -fabs(a)); }

template<class S>
inline static
void
SHFT(      S& a,
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



template<class S,class T, class Func>
void 
mnbrak(S&   ax,
       S&   bx,
       S&   cx,
       T&   fa,
       T&   fb,
       T&   fc,
       Func func)
  // given function func, and given distinct initial points ax and bx, this
  // routine searches in the downhill direction *defined by the function as
  // evaluated at the initial points) and returns new points ax, bx, cx which
  // bracket a minimum of the function. Also returned are the function values 
  // at the three points fa, fb and fc.
{	
  static const double GOLD   = 1.618034;
  static const double GLIMIT = 100.0;
  static const double TINY   = 1.0e-20;

 fa=func(ax);
 fb=func(bx);
 if (fb > fa) {
   CDS::swap(ax,bx);
   CDS::swap(fa,fb);
 }
 cx = bx + GOLD*(bx-ax);
 // cout << "nbrack: initial step: " << cx << ' ';
 fc = func(cx);
 while (fb > fc) {
   S r = (bx-ax) * (fb-fc);
   S q = (bx-cx) * (fb-fa);
   S u = bx - ((bx-cx)*q - (bx-ax)*r) /
	 (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
   S ulim = bx + GLIMIT*(cx-bx);
   T fu;
   if ( (bx-u)*(u-cx) > 0.0 ) {
     //     cout << "nbrack: parabolic step: " << u << ' ';
     fu = func(u);
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
     //     cout << "nbrack: min eval: " << u << ' ';
     fu = func(u);
   } else if ( (cx-u)*(u-ulim) > 0.0 ) {
     //     cout << "nbrack: hunt1 eval: " << u << ' ';
     T fu=func(u);
     if (fu < fc) {
       SHFT(bx,cx,u,cx+GOLD*(cx-bx));
       //       cout << "nbrack: hunt2 eval: " << u << ' ';
       SHFT(fb,fc,fu,func(u));
     }
   } else if ( (u-ulim)*(ulim-cx) >= 0.0 ) {
     u = ulim;
     //     cout << "nbrack: hunt3 eval: " << u << ' ';
     fu = func(u);
   } else {
     u  = cx + GOLD*(cx-bx);
     //     cout << "nbrack: hunt4 eval: " << u << ' ';
     fu = func(u);
   }
   SHFT(ax,bx,cx,u);
   SHFT(fa,fb,fc,fu);
 }
}

#endif /* __mnbrak_hh__ */
