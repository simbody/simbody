#ifndef __brent_hh__
#define __brent_hh__

#include <iostream>

template<class S,class Func>
S
brent(const S&      ax,
      const S&      bx,
      const S&      cx,
        Func    f,
        S&      xmin,
      const double& tol=1e-8,
      const int     ITMAX=100,
      const double& ZEPS=1e-10);


template<class S,class Func>
S
brent(const S&      ax,
      const S&      bx,
      const S&      cx,
        Func    f,
        S&      xmin,
      const double& tol,
      const int     ITMAX,
      const double& ZEPS)
{
    
 static const double CGOLD =  0.3819660;

 S a=((ax < cx) ? ax : cx);  //a and b bracket minimum w/
 S b=((ax > cx) ? ax : cx);  //  a<b
 S x = bx;
 S w = bx;
 S v = bx;
 S fw = f(x);
 S fv = fw;
 S fx = fw;

 S e=0;
 S d=0;
 for (int iter=1 ; iter<=ITMAX ; iter++) {
   S u,fu;
   S xm=0.5*(a+b);                            //midpoint of a and b
   S tol1 = tol*fabs(x)+ZEPS;
   S tol2= 2 * tol1;
   if (fabs(x-xm) <= (tol2-0.5*(b-a))) {      //stop if points close together
     xmin=x;
     return fx;
   }
   if (fabs(e) > tol1) {
     S r=(x-w)*(fx-fv);
     S q=(x-v)*(fx-fw);
     S p=(x-v)*q-(x-w)*r;
     q=2*(q-r);
     if (q > 0.0) p = -p;
     q=fabs(q);
     S etemp=e;
     e=d;
     if (fabs(p) >= fabs(0.5*q*etemp) || 
     p <= q*(a-x) || 
     p >= q*(b-x)                  ) {
       //       cout << "brent: taking parabolic step: ";
       d=CGOLD*(e=(x >= xm ? a-x : b-x));
     } else {
       //       cout << "brent: taking section step: ";
       d=p/q;
       u=x+d;
       if (u-a < tol2 || b-u < tol2)
     d=SIGN(tol1,xm-x);
     }
   } else {
     //     cout << "brent: taking section step: ";
     d=CGOLD*(e=(x >= xm ? a-x : b-x));
   }
   u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
   fu=f(u);
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
 std::cerr << "Too many iterations in BRENT\n";
 xmin=x;
 return fx;
}

#endif /* __brent_hh__ */
