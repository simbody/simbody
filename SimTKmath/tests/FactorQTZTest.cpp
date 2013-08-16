/* -------------------------------------------------------------------------- *
 *                          Simbody(tm): SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
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
        printf("\n  Estimated rank with default rcond  %d \n\n",qtz.getRank() );
        ASSERT( qtz.getRank() == 5 );

        qtz.factor(a, 0.01);
        qtz.solve( b, x );  // solve for x given a right hand side 
  
        printf("\n  Estimated rank with rcond = 0.01 : %d \n\n",qtz.getRank() );


        cout << " Overdetermined Double SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
        ASSERT((x-x_right).norm() < 0.001);

        FactorQTZ qtzCopy( qtz );
        Vector xc; // should get sized automatically to 5 by solve()
        qtzCopy.solve(b, xc );
        cout << " copy constructor SOLUTION:      " << xc << "  errnorm=" << (xc-x_right).norm() << endl;

        FactorQTZ qtzAssign = qtz;
        Vector xa; // should get sized automatically to 5 by solve()
        qtzCopy.solve(b, xa );
        cout << " copy assign SOLUTION:           " << xa << "  errnorm=" << (xa-x_right).norm() << endl;

        Matrix_<float> af(6,5); for (int i=0; i<6; ++i) for (int j=0; j<5; ++j) af(i,j)=(float)a(i,j);
        Vector_<float> bf(6); for (int i=0; i<6; ++i) bf[i] = (float)b[i];
        Vector_<float> xf_right(5); for (int i=0; i<5; ++i) xf_right[i] = (float)x_right[i];
        Vector_<float> xf; // should get sized automatically to 5 by solve()

        qtz.factor(af, (float)0.01);
        qtz.solve(bf,xf);

        cout << " Overdetermined Float SOLUTION:  " << xf << "  errnorm=" << (xf-xf_right).norm() << endl;
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

        Matrix bu2(3,2);
        bu2(0) = bu;
        bu2(1) = 2*bu;
        Matrix xu2; // should get sized 4x2 by solve

        FactorQTZ qtzu(au);  // perform QTZ factorization

        qtzu.solve( bu, xu );  // solve for x given a right hand side

        cout << " Underdetermined Double SOLUTION: " << xu << "  errnorm=" << (xu-xu_right).norm() << endl;
  
        qtzu.solve( bu2, xu2 );
        cout << " multiple rhs solution, double " << xu2 << endl;

        Matrix_<float> afu(3,4); for (int i=0; i<3; ++i) for (int j=0; j<4; ++j) afu(i,j)=(float)au(i,j);
        Vector_<float> bfu(3); for (int i=0; i<3; ++i) bfu[i] = (float)bu[i];
        Vector_<float> xfu_right(4); for (int i=0; i<4; ++i) xfu_right[i] = (float)xu_right[i];
        Vector_<float> xfu; // should get sized automatically to 4 by solve()

        Matrix_<float> bfu2(3,2);
        bfu2(0) = bfu;
        bfu2(1) = 2*bfu;
        Matrix_<float> xfu2; // should get sized 4x2 by solve

        FactorQTZ qtzfu(afu);  // perform QTZ factorization

        qtzfu.solve( bfu, xfu );  // solve for x given a right hand side
 
        cout << " Underdetermined Float SOLUTION: " << xfu << "  errnorm=" << (xfu-xfu_right).norm() << endl;

        qtzfu.solve( bfu2, xfu2 );
        cout << " multiple rhs solution, float " << xfu2 << endl;


       Real C[4] = { 1.0,   2.0,
              1.0,   3.0  };
        Matrix c(2,2, C);
        FactorQTZ cqtz(c);
        Matrix invQTZ;
        cqtz.inverse(invQTZ);
        cout << " FactorQTZ.inverse : " << endl;
        cout << invQTZ[0] << endl;
        cout << invQTZ[1] << endl;
  
        Real Z[4] = { 0.0,   0.0,
                     0.0,   0.0  };

        Matrix z(2,2, Z);
        FactorQTZ zqtz(z);
        Vector_<double> xz;
        Vector_<double> bz(2);
        bz(1) = bz(0) = 0.0;
        zqtz.solve( bz, xz );
        cout << " solve with mat all zeros : " << endl;
        for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");
        try {
            Matrix_<double> z0;
            FactorQTZ z0qtz(z0);
            Vector_<double> bz0(0);
            z0qtz.solve( bz0, xz );
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
    catch (...) {
        std::printf("FAILED: Unknown exception\n");
        return 1;
    }

    return 0;
}


