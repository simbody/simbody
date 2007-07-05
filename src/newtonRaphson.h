#ifndef SimTK_NEWTON_RAPHSON_H_
#define SimTK_NEWTON_RAPHSON_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
 * Contributors: Derived from IVM code written by Charles Schwieters.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"

#include <iostream>

class NewtonRaphson {
public:
    int    maxMin;
    int    maxIters;
    bool   verbose;
    Real zTol;  // minimization will quit if gradient norm is smaller
    std::ostream& errStream;
    String name;

    NewtonRaphson(const String& id, std::ostream& errStream=std::cerr)
      : maxMin( 20 ) , maxIters( 20 ), verbose(false), 
        zTol( 1e-8 ), errStream(errStream), name(id)
    { }

    // This will attempt to modify x to drive error b to below desiredTol, but will be
    // happy if it achieves only requiredTol (>= desiredTol). That is, as long as were
    // having an easy time of it we'll continue down to desiredTol, but we won't do
    // anything desperate once we're better than requiredTol.
    template<class CalcB,class CalcZ,class VecType>
    void calc(const Real& requiredTol, const Real& desiredTol, 
              VecType& x, CalcB calcB, CalcZ calcZ) const 
    {
        assert(requiredTol >= desiredTol);
        if ( verbose )
            errStream << "NewtonRaphson '" << name << "': start. Req="
            << requiredTol << " desired=" << desiredTol << "\n";
        VecType b = calcB(x);
        Real norm = std::sqrt(b.normSqr() / b.size());
        const Real zTolz = zTol*x.size();
        const Real zTol2 = zTolz*zTolz;
        Real onorm=norm;
        int  iters=0;
        if ( verbose )
            errStream << "NewtonRaphson '" << name << "': iter: " 
                      << iters << "  norm: " << norm << '\n';
        bool finished=(norm < desiredTol);
        while (!finished) {
            VecType ox = x;
            //std::cout << "NR: vars=" << x << std::endl; 
            //std::cout << "NR: errs=" << b << std::endl;
            
            VecType z = calcZ(b);
            x -= z;
            
            iters++;

            b = calcB(x);
            norm = std::sqrt(b.normSqr() / b.size());

            // If that made the norm worse, we're going to need to back off and feel 
            // our way around slowly here. But we won't bother if we've already met
            // requiredTol.
            if (norm > onorm) {
                if (onorm <= requiredTol) {
                    if (verbose) {
                        errStream << "NewtonRaphson '" << name 
                          << "': got required but not desired, giving up rather than gradient search\n";
                    }
                    x = ox; norm = onorm;   // back up to last good
                    finished = true;
                    continue;
                }

                if (verbose)
                    errStream << "NewtonRaphson '" << name << "': newton-Raphson failed."
                              << " Trying gradient search.\n";

                // Now we'll try a gradient search using up to maxMin iterations, halving
                // the step along the gradient each time. We'll stop at the first norm
                // which is an improvement over onorm.

                int mincnt = maxMin;
                do {
                    if ( mincnt < 1 ) 
                        SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                            "Too many gradient search steps taken");
                    z *= 0.5;
                    if ( z.normSqr() < zTol2 )  // TODO should be weighted norm of z
                        SimTK_THROW1(Exception::NewtonRaphsonFailure, "step along gradient too small");
                    x = ox - z;
                    b = calcB(x);
                    norm = std::sqrt(b.normSqr() / b.size());
                    if ( verbose )
                        errStream << "NewtonRaphson '" << name << "': gradient search in iter: " 
                                  << iters << "  norm: " << norm << '\n';
                    mincnt--;
                } while (norm >= onorm);
            }

            // At this point we know that norm is better than onorm.

            onorm = norm;

            if ( verbose )
                errStream << "NewtonRaphson '" << name << "': iter: " 
                          << iters << "  norm: " << norm << '\n';
            
            if (norm <= desiredTol)
                finished = true;
            else if (iters >= maxIters) {
                if (norm > requiredTol)
                    SimTK_THROW1(Exception::NewtonRaphsonFailure, "maxIters exceeded");
                if (verbose) {
                    errStream << "NewtonRaphson '" << name 
                      << "': got required but not desired, giving up after iter "
                      << iters << "\n";
                }
                finished = true;
            }
        }
    }
};

#endif // SimTK_NEWTON_RAPHSON_H_
