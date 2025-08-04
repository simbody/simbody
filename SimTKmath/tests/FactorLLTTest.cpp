/* -------------------------------------------------------------------------- *
 *                          Simbody(tm): SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-25 Stanford University and the Authors.        *
 * Authors: Alexander Beattie                                                 *
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
 * This is a test program which uses the FactorLLT  class to do an LLT 
 * factorization on a system of linear equations and then use the 
 * factored LLT matrix to solve for a particular right hand side 
 */

/*
The data for this test is from the following example:
URL:https://en.wikipedia.org/wiki/Cholesky_decomposition

Matrix:

     |   4   12  -16 |
A =  |  12   37  -43 |
     | -16  -43   98 |

Expected lower-triangular L:

     |  2   0   0 |
L =  |  6   1   0 |
     | -8   5   3 |

So that A = L * L^T

*/

#include "SimTKmath.h"

#include <iostream>
using std::cout; using std::endl;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

//  3x3 symmetric positive definite matrix
Real Adata[9] = {
    4,  12, -16,
    12,  37, -43,
    -16, -43,  98
};

Real Bdata[3] = {1, 2, 3};

// Checked with matlab
/*
format long g

A = [  4,  12, -16;
      12,  37, -43;
     -16, -43,  98 ];

b = [1; 2; 3];

% Perform Cholesky factorization and solve
R = chol(A); 
disp(R') % For L
norm(R'*R - A)
x = R\(R'\b)  
*/

Real Ldata[9] = {
    2, 0, 0,
    6, 1, 0,
    -8, 5, 3
};
Real Xdata[3] = { 28.5833333333333, -7.66666666666667,  1.33333333333333 };

int main () {
    try {

        Matrix A(3, 3, Adata);

        Vector b(3, Bdata);
        Vector x_right(3, Xdata);
        Vector x;

        // Perform LLT factorization and solve
        FactorLLT LLT(A);
        LLT.solve(b, x);
        cout << " Real SOLUTION: " << x << "  errnorm=" << (x-x_right).norm() << endl;
        ASSERT((x-x_right).norm() < 10*SignificantReal);

        Matrix L;
        LLT.getL(L);

        // Check that L is correct
        Matrix L_expected(3, 3, Ldata);
        Real L_error = (L - L_expected).norm();
        cout << " L matrix errnorm: " << L_error << endl;
        ASSERT(L_error < 10 * SignificantReal);

        // Reconstruct A from L * L^T and check accuracy
        Matrix A_reconstructed = L * ~L;
        Real errNorm = (A_reconstructed - A).norm();
        cout << " Reconstructed A = L * L^T errnorm=" << errNorm << endl;
        ASSERT(errNorm < 10 * SignificantReal);

    } catch (const std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }

    return 0;
}


