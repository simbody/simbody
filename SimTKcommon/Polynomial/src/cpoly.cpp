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
 * 01.01    HVE/20021101      Initial release
 * 01.02    PE/20070808       Converted to a class, templatized  
 * 01.03    MAS/20130410      Minor changes to match RPoly changes  
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

#include <cmath>
#include <limits>
#include "cpoly.h"

namespace SimTK {

/* define version string */
template<class T>
char CPoly<T>::_V_[] = "@(#)cpoly.h 01.03 -- Copyright (C) Henrik Vestermark";

template<class T>
int CPoly<T>::findRoots( const T *opr, const T *opi, int degree, T *zeror, T *zeroi ) 
   {
   int cnt1, cnt2, idnn2, i, conv;
   T xx, yy, cosr, sinr, smalno, base, xxx, zr, zi, bnd;

   mcon( &eta, &infin, &smalno, &base );
   are = eta;
   mre = (T) (2.0 * sqrt( 2.0 ) * eta);
   xx = (T) 0.70710678;
   yy = -xx;
   cosr = (T) -0.060756474;
   sinr = (T) -0.99756405;
   nn = degree;  

   // Algorithm fails if the leading coefficient is zero, or degree is zero.
   if( nn < 1 || (opr[ 0 ] == 0 && opi[ 0 ] == 0) )
      return -1;


   // Remove the zeros at the origin if any
   while( opr[ nn ] == 0 && opi[ nn ] == 0 )
      {
      idnn2 = degree - nn;
      zeror[ idnn2 ] = 0;
      zeroi[ idnn2 ] = 0;
      nn--;
      }

   // sherm 20130410: If all coefficients but the leading one were zero, then
   // all solutions are zero; should be a successful (if boring) return.
   if (nn == 0) 
      return degree;

   // Allocate arrays
   pr = new T [ degree+1 ];
   pi = new T [ degree+1 ];
   hr = new T [ degree+1 ];
   hi = new T [ degree+1 ];
   qpr= new T [ degree+1 ];
   qpi= new T [ degree+1 ];
   qhr= new T [ degree+1 ];
   qhi= new T [ degree+1 ];
   shr= new T [ degree+1 ];
   shi= new T [ degree+1 ];

   // Make a copy of the coefficients
   for( i = 0; i <= nn; i++ )
      {
      pr[ i ] = opr[ i ];
      pi[ i ] = opi[ i ];
      shr[ i ] = cmod( pr[ i ], pi[ i ] );
      }

   // Scale the polynomial
   bnd = scale( nn, shr, eta, infin, smalno, base );
   if( bnd != 1 )
      for( i = 0; i <= nn; i++ )
         {
         pr[ i ] *= bnd;
         pi[ i ] *= bnd;
         }

search: 
   if( nn <= 1 )
      {
      cdivid( -pr[ 1 ], -pi[ 1 ], pr[ 0 ], pi[ 0 ], &zeror[ degree-1 ], &zeroi[ degree-1 ] );
      goto finish;
      }

   // Calculate bnd, alower bound on the modulus of the zeros
   for( i = 0; i<= nn; i++ )
      shr[ i ] = cmod( pr[ i ], pi[ i ] );

   cauchy( nn, shr, shi, &bnd );
   
   // Outer loop to control 2 Major passes with different sequences of shifts
   for( cnt1 = 1; cnt1 <= 2; cnt1++ )
      {
      // First stage  calculation , no shift
      noshft( 5 );

      // Inner loop to select a shift
      for( cnt2 = 1; cnt2 <= 9; cnt2++ )
         {
         // Shift is chosen with modulus bnd and amplitude rotated by 94 degree from the previous shif
         xxx = cosr * xx - sinr * yy;
         yy = sinr * xx + cosr * yy;
         xx = xxx;
         sr = bnd * xx;
         si = bnd * yy;

         // Second stage calculation, fixed shift
         fxshft( 10 * cnt2, &zr, &zi, &conv );
         if( conv )
            {
            // The second stage jumps directly to the third stage ieration
            // If successful the zero is stored and the polynomial deflated
            idnn2 = degree - nn;
            zeror[ idnn2 ] = zr;
            zeroi[ idnn2 ] = zi;
            nn--;
            for( i = 0; i <= nn; i++ )
               {
               pr[ i ] = qpr[ i ];
               pi[ i ] = qpi[ i ];
               }
            goto search;
            }
         // If the iteration is unsuccessful another shift is chosen
         }
      // if 9 shifts fail, the outer loop is repeated with another sequence of shifts
      }

   // The zerofinder has failed on two major passes
   // return empty handed with the number of roots found (less than the original degree)
   degree -= nn;

finish:
   // Deallocate arrays
   delete [] pr;
   delete [] pi;
   delete [] hr;
   delete [] hi;
   delete [] qpr;
   delete [] qpi;
   delete [] qhr;
   delete [] qhi;
   delete [] shr;
   delete [] shi;

   return degree;       
   }


// COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
// POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
//
template<class T>
void CPoly<T>::noshft( const int l1 )
   {
   int i, j, jj, n, nm1;
   T xni, t1, t2;

   n = nn;
   nm1 = n - 1;
   for( i = 0; i < n; i++ )
      {
      xni = (T) (nn - i);
      hr[ i ] = xni * pr[ i ] / n;
      hi[ i ] = xni * pi[ i ] / n;
      }
   for( jj = 1; jj <= l1; jj++ )
      {
      if( cmod( hr[ n - 1 ], hi[ n - 1 ] ) > eta * 10 * cmod( pr[ n - 1 ], pi[ n - 1 ] ) )
         {
         cdivid( -pr[ nn ], -pi[ nn ], hr[ n - 1 ], hi[ n - 1 ], &tr, &ti );
         for( i = 0; i < nm1; i++ )
            {
            j = nn - i - 1;
            t1 = hr[ j - 1 ];
            t2 = hi[ j - 1 ];
            hr[ j ] = tr * t1 - ti * t2 + pr[ j ];
            hi[ j ] = tr * t2 + ti * t1 + pi[ j ];
            }
         hr[ 0 ] = pr[ 0 ];
         hi[ 0 ] = pi[ 0 ];
         }
      else
         {
         // If the constant term is essentially zero, shift H coefficients
         for( i = 0; i < nm1; i++ )
            {
            j = nn - i - 1;
            hr[ j ] = hr[ j - 1 ];
            hi[ j ] = hi[ j - 1 ];
            }
         hr[ 0 ] = 0;
         hi[ 0 ] = 0;
         }
      }
   }

// COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
// INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
// APPROXIMATE ZERO IF SUCCESSFUL.
// L2 - LIMIT OF FIXED SHIFT STEPS
// ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
// CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
//
template<class T>
void CPoly<T>::fxshft( const int l2, T *zr, T *zi, int *conv )
   {
   int i, j, n;
   int test, pasd, bol;
   T otr, oti, svsr, svsi;

   n = nn;
   polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
   test = 1;
   pasd = 0;

   // Calculate first T = -P(S)/H(S)
   calct( &bol );

   // Main loop for second stage
   for( j = 1; j <= l2; j++ )
      {
      otr = tr;
      oti = ti;

      // Compute the next H Polynomial and new t
      nexth( bol );
      calct( &bol );
      *zr = sr + tr;
      *zi = si + ti;

      // Test for convergence unless stage 3 has failed once or this
      // is the last H Polynomial
      if( !( bol || !test || j == 12 ) )
         {
         if( cmod( tr - otr, ti - oti ) < 0.5 * cmod( *zr, *zi ) )
            {
            if( pasd )
               {
               // The weak convergence test has been passed twice, start the third stage
               // Iteration, after saving the current H polynomial and shift
               for( i = 0; i < n; i++ )
                  {
                  shr[ i ] = hr[ i ];
                  shi[ i ] = hi[ i ];
                  }
               svsr = sr;
               svsi = si;
               vrshft( 10, zr, zi, conv );
               if( *conv ) return;

               //The iteration failed to converge. Turn off testing and restore h,s,pv and T
               test = 0;
               for( i = 0; i < n; i++ )
                  {
                  hr[ i ] = shr[ i ];
                  hi[ i ] = shi[ i ];
                  }
               sr = svsr;
               si = svsi;
               polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
               calct( &bol );
               continue;
               }
            pasd = 1;
            }
         else
            pasd = 0;
         }
      }

   // Attempt an iteration with final H polynomial from second stage
   vrshft( 10, zr, zi, conv );
   }

// CARRIES OUT THE THIRD STAGE ITERATION.
// L3 - LIMIT OF STEPS IN STAGE 3.
// ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
//           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
// CONV    -  .TRUE. IF ITERATION CONVERGES
//
template<class T>
void CPoly<T>::vrshft( const int l3, T *zr, T *zi, int *conv )
   {
   int b, bol;
   int i, j;
   T mp, ms, omp, relstp, r1, r2, tp;

   *conv = 0;
   b = 0;
   sr = *zr;
   si = *zi;

   // Main loop for stage three
   for( i = 1; i <= l3; i++ )
      {
      // Evaluate P at S and test for convergence
      polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
      mp = cmod( pvr, pvi );
      ms = cmod( sr, si );
      if( mp <= 20 * errev( nn, qpr, qpi, ms, mp, are, mre ) )
         {
         // Polynomial value is smaller in value than a bound onthe error
         // in evaluationg P, terminate the ietartion
         *conv = 1;
         *zr = sr;
         *zi = si;
         return;
         }
      if( i != 1 )
         {
         if( !( b || mp < omp || relstp >= 0.05 ) )
            {
            // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed 
            // shift steps into the cluster to force one zero to dominate
            tp = relstp;
            b = 1;
            if( relstp < eta ) tp = eta;
            r1 = sqrt( tp );
            r2 = sr * ( 1 + r1 ) - si * r1;
            si = sr * r1 + si * ( 1 + r1 );
            sr = r2;
            polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
            for( j = 1; j <= 5; j++ )
               {
               calct( &bol );
               nexth( bol );
               }
            omp = infin;
            goto _20;
            }
         
         // Exit if polynomial value increase significantly
         if( mp *0.1 > omp ) return;
         }

      omp = mp;

      // Calculate next iterate
_20:  calct( &bol );
      nexth( bol );
      calct( &bol );
      if( !bol )
         {
         relstp = cmod( tr, ti ) / cmod( sr, si );
         sr += tr;
         si += ti;
         }
      }
   }

// COMPUTES  T = -P(S)/H(S).
// BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
template<class T>
void CPoly<T>::calct( int *bol )
   {
   int n;
   T hvr, hvi;

   n = nn;

   // evaluate h(s)
   polyev( n - 1, sr, si, hr, hi, qhr, qhi, &hvr, &hvi );
   *bol = cmod( hvr, hvi ) <= are * 10 * cmod( hr[ n - 1 ], hi[ n - 1 ] ) ? 1 : 0;
   if( !*bol )
      {
      cdivid( -pvr, -pvi, hvr, hvi, &tr, &ti );
      return;
      }

   tr = 0;
   ti = 0;
   }

// CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
// BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
//
template<class T>
void CPoly<T>::nexth( const int bol )
   {
   int j, n;
   T t1, t2;

   n = nn;
   if( !bol )
      {
      for( j = 1; j < n; j++ )
         {
         t1 = qhr[ j - 1 ];
         t2 = qhi[ j - 1 ];
         hr[ j ] = tr * t1 - ti * t2 + qpr[ j ];
         hi[ j ] = tr * t2 + ti * t1 + qpi[ j ];
         }
      hr[ 0 ] = qpr[ 0 ];
      hi[ 0 ] = qpi[ 0 ];
      return;
      }

   // If h[s] is zero replace H with qh
   for( j = 1; j < n; j++ )
      {
      hr[ j ] = qhr[ j - 1 ];
      hi[ j ] = qhi[ j - 1 ];
      }
   hr[ 0 ] = 0;
   hi[ 0 ] = 0;
   }

// EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
// PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
//  
template<class T>
void CPoly<T>::polyev( const int nn, const T sr, const T si, const T pr[], const T pi[], T qr[], T qi[], T *pvr, T *pvi )  
   {
   int i;
   T t;

   qr[ 0 ] = pr[ 0 ];
   qi[ 0 ] = pi[ 0 ];
   *pvr = qr[ 0 ];
   *pvi = qi[ 0 ];

   for( i = 1; i <= nn; i++ )
      {
      t = ( *pvr ) * sr - ( *pvi ) * si + pr[ i ];
      *pvi = ( *pvr ) * si + ( *pvi ) * sr + pi[ i ];
      *pvr = t;
      qr[ i ] = *pvr;
      qi[ i ] = *pvi;
      }
   }

// BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER RECURRENCE.
// QR,QI - THE PARTIAL SUMS
// MS    -MODULUS OF THE POINT
// MP    -MODULUS OF POLYNOMIAL VALUE
// ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION
//
template<class T>
T CPoly<T>::errev( const int nn, const T qr[], const T qi[], const T ms, const T mp, const T are, const T mre )
   {
   int i;
   T e;

   e = cmod( qr[ 0 ], qi[ 0 ] ) * mre / ( are + mre );
   for( i = 0; i <= nn; i++ )
      e = e * ms + cmod( qr[ i ], qi[ i ] );

   return e * ( are + mre ) - mp * mre;
   }

// CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
// POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
//
template<class T>
void CPoly<T>::cauchy( const int nn, T pt[], T q[], T *fn_val )
   {
   int i, n;
   T x, xm, f, dx, df;

   pt[ nn ] = -pt[ nn ];

   // Compute upper estimate bound
   n = nn;
   x = exp( log( -pt[ nn ] ) - log( pt[ 0 ] ) ) / n;
   if( pt[ n - 1 ] != 0 )
      {
      // Newton step at the origin is better, use it
      xm = -pt[ nn ] / pt[ n - 1 ];
      if( xm < x ) x = xm;
      }

   // Chop the interval (0,x) until f < 0
   while(1)
      {
      xm = x * (T) 0.1;
      f = pt[ 0 ];
      for( i = 1; i <= nn; i++ )
         f = f * xm + pt[ i ];
      if( f <= 0 )
         break;
      x = xm;
      }
   dx = x;
   
   // Do Newton iteration until x converges to two decimal places
   while( std::abs( dx / x ) > 0.005 )
      {
      q[ 0 ] = pt[ 0 ];
      for( i = 1; i <= nn; i++ )
         q[ i ] = q[ i - 1 ] * x + pt[ i ];
      f = q[ nn ];
      df = q[ 0 ];
      for( i = 1; i < n; i++ )
         df = df * x + q[ i ];
      dx = f / df;
      x -= dx;
      }

   *fn_val = x;
   }

// RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.
// THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
// INTERFERING WITH THE CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE
// BASE.
// PT - MODULUS OF COEFFICIENTS OF P
// ETA, INFIN, SMALNO, BASE - CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC.
//
template<class T>
T CPoly<T>::scale( const int nn, const T pt[], const T eta, const T infin, const T smalno, const T base )
   {
   int i, l;
   T hi, lo, max, min, x, sc;
   T fn_val;

   // Find largest and smallest moduli of coefficients
   hi = sqrt( infin );
   lo = smalno / eta;
   max = 0;
   min = infin;

   for( i = 0; i <= nn; i++ )
      {
      x = pt[ i ];
      if( x > max ) max = x;
      if( x != 0 && x < min ) min = x;
      }

   // Scale only if there are very large or very small components
   fn_val = 1;
   if( min >= lo && max <= hi ) return fn_val;
   x = lo / min;
   if( x <= 1 )
      sc = 1 / ( sqrt( max )* sqrt( min ) );
   else
      {
      sc = x;
      if( infin / sc > max ) sc = 1;
      }
   l = (int)( log( sc ) / log(base ) + 0.5 );
   fn_val = pow( base , l );
   return fn_val;
   }

// COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
//
template<class T>
void CPoly<T>::cdivid( const T ar, const T ai, const T br, const T bi, T *cr, T *ci )
   {
   T r, d, t, infin;

   if( br == 0 && bi == 0 )
      {
      // Division by zero, c = infinity
      mcon( &t, &infin, &t, &t );
      *cr = infin;
      *ci = infin;
      return;
      }

   if( std::abs( br ) < std::abs( bi ) )
      {
      r = br/ bi;
      d = bi + r * br;
      *cr = ( ar * r + ai ) / d;
      *ci = ( ai * r - ar ) / d;
      return;
      }

   r = bi / br;
   d = br + r * bi;
   *cr = ( ar + ai * r ) / d;
   *ci = ( ai - ar * r ) / d;
   }

// MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
//
template<class T>
T CPoly<T>::cmod( const T r, const T i )
   {
   T ar, ai;

   ar = std::abs( r );
   ai = std::abs( i );
   if( ar < ai )
      return ai * sqrt( (T) 1.0 + pow( ( ar / ai ), (T) 2.0 ) );

   if( ar > ai )
      return ar * sqrt( (T) 1.0 + pow( ( ai / ar ),(T) 2.0 ) );

   return ar * sqrt( (T) 2.0 );
   }

// MCON PROVIDES MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE PROGRAM.
// THE USER MAY EITHER SET THEM DIRECTLY OR USE THE STATEMENTS BELOW TO
// COMPUTE THEM. THE MEANING OF THE FOUR CONSTANTS ARE -
// ETA       THE MAXIMUM RELATIVE REPRESENTATION ERROR WHICH CAN BE DESCRIBED
//           AS THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
//           1.0_dp + ETA &gt; 1.0.
// INFINY    THE LARGEST FLOATING-POINT NUMBER
// SMALNO    THE SMALLEST POSITIVE FLOATING-POINT NUMBER
// BASE      THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED
//
template<class T>
void CPoly<T>::mcon( T *eta, T *infiny, T *smalno, T *base )
   {
   *base = FLT_RADIX;
   *eta = std::numeric_limits<T>::epsilon();
   *infiny = std::numeric_limits<T>::max();
   *smalno = std::numeric_limits<T>::min();
   }

template class CPoly<float>;
template class CPoly<double>;

} // namespace SimTK
