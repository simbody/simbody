#ifndef SimTK_SimTKCOMMON_RPOLY_H_
#define SimTK_SimTKCOMMON_RPOLY_H_

/*      rpoly.cpp -- Jenkins-Traub real polynomial root finder.
 *
 *      Written by C. Bond, with minor changes by Peter Eastman.
 *      This file is in the public domain.
 *
 *      Translation of TOMS493 from FORTRAN to C. This
 *      implementation of Jenkins-Traub partially adapts
 *      the original code to a C environment by restruction
 *      many of the 'goto' controls to better fit a block
 *      structured form. It also eliminates the global memory
 *      allocation in favor of local, dynamic memory management.
 *
 *      The calling conventions are slightly modified to return
 *      the number of roots found as the function value.
 *
 *      INPUT:
 *      op - vector of coefficients in order of
 *              decreasing powers.
 *      degree - integer degree of polynomial
 *
 *      OUTPUT:
 *      zeror,zeroi - output vectors of the
 *              real and imaginary parts of the zeros.
 *
 *      RETURN:
 *      returnval:   -1 if leading coefficient is zero, otherwise
 *                  number of roots found. 
 */

namespace SimTK {

template <class T>
class RPoly {
public:
    int findRoots(T *op, int degree, T *zeror, T *zeroi);
private:
    void quad(T a,T b1,T c,T *sr,T *si,
            T *lr,T *li);
    void fxshfr(int l2, int *nz);
    void quadit(T *uu,T *vv,int *nz);
    void realit(T *sss, int *nz, int *iflag);
    void calcsc(int *type);
    void nextk(int *type);
    void newest(int type,T *uu,T *vv);
    void quadsd(int n,T *u,T *v,T *p,T *q,
            T *a,T *b);
    T *p,*qp,*k,*qk,*svk;
    T sr,si,u,v,a,b,c,d,a1,a2;
    T a3,a6,a7,e,f,g,h,szr,szi,lzr,lzi;
    T eta,are,mre;
    int n,nn,nmi,zerok;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_RPOLY_H_
