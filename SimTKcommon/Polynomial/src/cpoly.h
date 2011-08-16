#ifndef SimTK_SimTKCOMMON_CPOLY_H_
#define SimTK_SimTKCOMMON_CPOLY_H_

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002
 *                       Henrik Vestermark
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   Permission to use, copy, distribute, and sell this software and its
 *   documentation for any purpose is hereby granted without fee, provided:
 *   THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY WARRANTY
 *   OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL
 *   Henrik Vestermark or Future Team Aps, BE LIABLE FOR ANY SPECIAL,
 *   INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR ANY DAMAGES
 *   WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER OR NOT
 *   ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF LIABILITY,
 *   ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 *   SOFTWARE.
 *
 *******************************************************************************
 *
 *
 * Module name     :   cpoly.cpp
 * Module ID Nbr   :   
 * Description     :   cpoly.cpp -- Jenkins-Traub real polynomial root finder.
 *                     Translation of TOMS493 from FORTRAN to C. This
 *                     implementation of Jenkins-Traub partially adapts
 *                     the original code to a C environment by restruction
 *                     many of the 'goto' controls to better fit a block
 *                     structured form. It also eliminates the global memory
 *                     allocation in favor of local, dynamic memory management.
 *
 *                     The calling conventions are slightly modified to return
 *                     the number of roots found as the function value.
 *
 *                     INPUT:
 *                     opr - vector of real coefficients in order of
 *                          decreasing powers.
 *                     opi - vector of imaginary coefficients in order of
 *                          decreasing powers.
 *                     degree - integer degree of polynomial
 *
 *                     OUTPUT:
 *                     zeror,zeroi - output vectors of the
 *                          real and imaginary parts of the zeros.
 *                            to be consistent with rpoly.cpp the zeros is inthe index
 *                            [0..max_degree-1]
 *
 *                     RETURN:
 *                     returnval:   -1 if leading coefficient is zero, otherwise
 *                          number of roots found. 
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version  Author/Date     Description of changes
 * -------  -----------     ----------------------
 * 01.01    HVE/021101      Initial release
 * 01.02    PE/8Aug07       Converted to a class, templatized  
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

#include <float.h>

namespace SimTK {

template<class T>
class CPoly {
public:
    static char _V_[];
    int findRoots( const T *opr, const T *opi, int degree, T *zeror, T *zeroi );
private:
    T sr, si, tr, ti, pvr, pvi, are, mre, eta, infin;
    int nn;
    T *pr, *pi, *hr, *hi, *qpr, *qpi, *qhr, *qhi, *shr, *shi; 
    void noshft( const int l1 );
    void fxshft( const int l2, T *zr, T *zi, int *conv );
    void vrshft( const int l3, T *zr, T *zi, int *conv );
    void calct( int *bol );
    void nexth( const int bol );
    void polyev( const int nn, const T sr, const T si, const T pr[], const T pi[], T qr[], T qi[], T *pvr, T *pvi );
    T errev( const int nn, const T qr[], const T qi[], const T ms, const T mp, const T are, const T mre );
    void cauchy( const int nn, T pt[], T q[], T *fn_val );
    T scale( const int nn, const T pt[], const T eta, const T infin, const T smalno, const T base );
    void cdivid( const T ar, const T ai, const T br, const T bi, T *cr, T *ci );
    T cmod( const T r, const T i );
    void mcon( T *eta, T *infiny, T *smalno, T *base );
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_CPOLY_H_
