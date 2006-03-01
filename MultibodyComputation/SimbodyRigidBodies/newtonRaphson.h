#ifndef NEWTON_RAPHSON_H_
#define NEWTON_RAPHSON_H_

#include "simbody/internal/SimbodyCommon.h"

#include <iostream>

class NewtonRaphson {
public:
    int    maxMin;
    int    maxIters;
    bool   verbose;
    double tol;
    double zTol;  // minimization will quit if gradient norm is smaller
    std::ostream& errStream;

    NewtonRaphson(std::ostream& errStream=std::cerr)
      : maxMin( 20 ) , maxIters( 20 ), verbose(0), tol(1e-8), 
        zTol( 1e-8 ), errStream(errStream)
    { }

    template<class CalcB,class CalcZ,class VecType>
    void calc(VecType& x, CalcB calcB, CalcZ calcZ) 
    {
        if ( verbose )
            errStream << "NewtonRaphson: start.\n";
        VecType b = calcB(x);
        double norm = b.norm() / x.size();
        const double zTolz = zTol*x.size();
        const double zTol2 = zTolz*zTolz;
        double onorm=norm;
        int  iters=0;
        if ( verbose )
            errStream << "NewtonRaphson: iter: " 
                      << iters << "  norm: " << norm << '\n';
        bool finished=(norm < tol);
        while (!finished) {
            VecType ox = x;
            
            VecType z = calcZ(b);
            x += z;
            
            iters++;

            b = calcB(x);
            norm = b.norm() / x.size();

            if ( norm > onorm && verbose)
                errStream << "NewtonRaphson: newton-Raphson failed."
                          << " Trying gradient search.\n";

            int mincnt = maxMin;
            z *= -1.0;
            while ( norm > onorm ) {
                if ( mincnt < 1 ) 
                    SIMTK_THROW1(Exception::NewtonRaphsonFailure, "too many minimization steps taken");
                z *= 0.5;
                if ( z.normSqr() < zTol2 )
                    SIMTK_THROW1(Exception::NewtonRaphsonFailure, "gradient too small");
                x = ox + z;
                b = calcB(x);
                norm = b.norm() / x.size();
                if ( verbose )
                    errStream << "NewtonRaphson: iter: " 
                              << iters << "  norm: " << norm << '\n';
                mincnt--;
            }

            onorm = norm;

            if ( verbose )
                errStream << "NewtonRaphson: iter: " 
                          << iters << "  norm: " << norm << '\n';
            
            if (norm < tol)
                finished=1;
            if (iters > maxIters) 
                SIMTK_THROW1(Exception::NewtonRaphsonFailure, "maxIters exceeded");
        }
    }
};

#endif // NEWTON_RAPHSON_H_
