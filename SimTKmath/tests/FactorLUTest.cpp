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


Solves:

Ax = B,

where A is the general matrix


      1.80   2.88   2.05   -0.89           9.52
A =   5.25  -2.95  -0.95   -3.80   and B = 24.35
      1.58  -2.69  -2.90   -1.04           0.77
     -1.11  -0.66  -0.59    0.80          -6.22


   Solution
 x =     1.0000    -1.0000     3.0000    -5.0000


 LU factorization:
             1          2          3          4
 1      5.2500    -2.9500    -0.9500    -3.8000
 2      0.3429     3.8914     2.3757     0.4129
 3      0.3010    -0.4631    -1.5139     0.2948
 4     -0.2114    -0.3299     0.0047     0.1314

 Pivot indices
             2          2          3          4

*/

#include "SimTKmath.h"

#include <iostream>
using std::cout; using std::endl;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

Real A[16] = { 1.80,   2.88,   2.05,   -0.89,
               5.25,  -2.95,  -0.95,   -3.80,
               1.58,  -2.69,  -2.90,   -1.04,
              -1.11,  -0.66,  -0.59,    0.80  };

Real B[4] =  { 9.52, 24.35,  0.77, -6.22 };

Real X[4] =  { 1.,   -1.,    3.,   -5.   };

int main () {
    try { 
            // Default precision (Real, normally double) test.
        Matrix a(4,4, A);
        Vector b(4, B);
        Vector x_right(4, X);
        Vector x; // should get sized automatically to 4 by solve()

        FactorLU lu(a);  // perform LU factorization 

        lu.solve( b, x );  // solve for x given a right hand side 

        cout << " Real SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
        ASSERT((x-x_right).norm() < 10*SignificantReal);

            // float test

        Matrix_<float> af(4,4); for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) af(i,j)=(float)a(i,j);
        Vector_<float> bf(4); for (int i=0; i<4; ++i) bf[i] = (float)b[i];
        Vector_<float> xf_right(4); for (int i=0; i<4; ++i) xf_right[i] = (float)x_right[i];
        Vector_<float> xf; // should get sized automatically to 4 by solve()

        FactorLU luf;
        luf.factor(af);
        luf.solve(bf, xf);

        cout << " float SOLUTION: " << xf << "  errnorm=" << (xf-xf_right).norm() << endl;
        const float SignificantFloat = NTraits<float>::getSignificant();
        ASSERT((xf-xf_right).norm() < 10*SignificantFloat);

        luf.factor(a);
        lu.solve( b, x );  // solve for x given a right hand side 
        cout << " Real SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
        ASSERT((x-x_right).norm() < 10*SignificantReal);
        
        Real C[4] = { 1.0,   2.0,
                      1.0,   3.0  };
        Matrix c(2,2, C);
        FactorLU clu(c);
        Matrix invC;
        clu.inverse(invC);
        cout << "Inverse c: " << endl;
        cout << invC[0] << endl;
        cout << invC[1] << endl;
        Real Z[4] = { 0.0,   0.0,
                     0.0,   0.0  };
        Matrix z(2,2, Z);
        FactorLU zlu(z);
        Vector_<double> xz;
        Vector_<double> bz(2);
        bz(1) = bz(0) = 0.0;
        zlu.solve( bz, xz );
        cout << " solve with mat all zeros : " << endl;
        for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");
   
        try {
            Matrix_<double> z0;
            FactorLU z0lu(z0);
            Vector_<double> bz0(0);
            z0lu.solve( bz0, xz );
            cout << " solve with mat(0,0) : " << endl;
            for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");
        } catch (const std::exception& e) {
             cout << "(EXPECTED EXCEPTION) NULL matrix test: " 
                 << e.what() << endl;
        }
    } 
    catch (const std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }

    return 0;
}


