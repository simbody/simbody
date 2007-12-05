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
 * This is a test program which uses the FactorLU  class to do an LU 
 * factorization on a system of linear equations and then use the 
 * factored LU matrix to solve for a particular right hand side 
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

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

  double A[30] = { -0.09,   0.14,  -0.46,    0.68,   1.29,       
                   -1.56,   0.20,   0.29,    1.09,   0.51,        
                   -1.48,  -0.43,   0.89,   -0.71,  -0.96,   
                   -1.09,   0.84,   0.77,    2.11,  -1.27,       
                    0.08,   0.55,  -1.13,    0.14,   1.74,        
                   -1.59,  -0.72,   1.06,    1.24,   0.34  };  

  double B[6] =  { 7.4, 4.2, -8.3, 1.8, 8.6, 2.1 };

int main () {
    try { 
        int i,j;
        Matrix_<double> a(6,5, A);
        Vector_<double> b(6, B);
        Vector_<double> x(6);

        FactorQTZ qtz(a);  // perform LU factorization 

        qtz.solve( b, x );  // solve for x given a right hand side 

        printf(" SOLUTION: "); for(i=0;i<5;i++) printf("%f ",x(i)); printf(" \n" );
        
        return 0;
    } 
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}


