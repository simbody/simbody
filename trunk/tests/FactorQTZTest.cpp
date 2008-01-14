/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Jack Middleton                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**@file
 * This is a test program which uses the FactorQTZ  class to do an QTZ 
 * factorization on a system of linear equations and then use the 
 * factored QTZ matrix to solve a find a least squares solution 
 * for a particular right hand side 
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves the least squares problem:

Ax = B,

where A is the general matrix


     -0.09   0.14  -0.46    0.68   1.29        7.4
     -1.56   0.20   0.29    1.09   0.51        4.2
A =  -1.48  -0.43   0.89   -0.71  -0.96    B= -8.3
     -1.09   0.84   0.77    2.11  -1.27        1.8
      0.08   0.55  -1.13    0.14   1.74        8.6
     -1.59  -0.72   1.06    1.24   0.34        2.1

The default tolerance of 0.01 is used to determine the effective rank of A


SOLUTION = 
0.6344     0.9699    -1.4402     3.3678     3.3992

estimated rank = 4

*/



//#define SimTK_USE_STATIC_LIBRARIES

#include "SimTKcommon.h"
#include "SimTKmath.h"

#include <cstdio>
#include <cassert>
#include <iostream>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

Real A[30] = { -0.09,   0.14,  -0.46,    0.68,   1.29,       
                   -1.56,   0.20,   0.29,    1.09,   0.51,        
                   -1.48,  -0.43,   0.89,   -0.71,  -0.96,   
                   -1.09,   0.84,   0.77,    2.11,  -1.27,       
                    0.08,   0.55,  -1.13,    0.14,   1.74,        
                   -1.59,  -0.72,   1.06,    1.24,   0.34  };  

Real B[6] =  { 7.4, 4.2, -8.3, 1.8, 8.6, 2.1 };
Real X[5] =  { 0.6344, 0.9699, -1.4402, 3.3678,  3.3992 };

int main () {
    try { 
           // Default precision (Real, normally double) test.

        Matrix a(6,5, A);
        Vector b(6, B);
        Vector x_right(5, X);
        Vector x; // should get sized automatically to 5 by solve()

        FactorQTZ qtz;  // perform QTZ factorization 

        qtz.factor(a);
        qtz.solve( b, x );  // solve for x given a right hand side 


        cout << " Overdetermined Double SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
//        ASSERT((x-x_right).norm() < 10*SignificantReal);
        ASSERT((x-x_right).norm() < 0.001);

        Matrix_<float> af(6,5); for (int i=0; i<6; ++i) for (int j=0; j<5; ++j) af(i,j)=(float)a(i,j);
        Vector_<float> bf(6); for (int i=0; i<6; ++i) bf[i] = (float)b[i];
        Vector_<float> xf_right(5); for (int i=0; i<5; ++i) xf_right[i] = (float)x_right[i];
        Vector_<float> xf; // should get sized automatically to 5 by solve()

          qtz.factor(af);
          qtz.solve(bf,xf);

        cout << " Overdetermined Float SOLUTION: " << xf << "  errnorm=" << (xf-xf_right).norm() << endl;
        const float SignificantFloat = NTraits<float>::getSignificant();
        ASSERT((xf-xf_right).norm() < 0.001);

        // Underdetermined case adapted from 
        // http://idlastro.gsfc.nasa.gov/idl_html_help/LA_LEAST_SQUARES.html
        
        Real Au[12] = { 2,     5,     3,     4,
                        7,     1,     3,     5,
                        4,     3,     6,     2   };
        Real Bu[3] = { 3,     1,     6 };
        Real Xu[4] = { -0.0376844,     0.350628,    0.986164,   -0.409066 };
        Matrix au(3, 4, Au);
        Vector bu(3, Bu);
        Vector xu_right(4, Xu);
        Vector xu; // should get sized automatically to 4 by solve()

        FactorQTZ qtzu(au);  // perform QTZ factorization

        qtzu.solve( bu, xu );  // solve for x given a right hand side

        cout << " Underdetermined Double SOLUTION: " << xu << "  errnorm=" << (xu-xu_right).norm() << endl;
  
        Matrix_<float> afu(3,4); for (int i=0; i<3; ++i) for (int j=0; j<4; ++j) afu(i,j)=(float)au(i,j);
        Vector_<float> bfu(3); for (int i=0; i<3; ++i) bfu[i] = (float)bu[i];
        Vector_<float> xfu_right(4); for (int i=0; i<4; ++i) xfu_right[i] = (float)xu_right[i];
        Vector_<float> xfu; // should get sized automatically to 4 by solve()

        FactorQTZ qtzfu(afu);  // perform QTZ factorization

        qtzfu.solve( bfu, xfu );  // solve for x given a right hand side
 
        cout << " Underdetermined Float SOLUTION: " << xfu << "  errnorm=" << (xfu-xfu_right).norm() << endl;
  
        return 0;
    } 
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}


