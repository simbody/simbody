#ifndef CDS_NEWTON_RAPHSON_H_
#define CDS_NEWTON_RAPHSON_H_

#include "cdsMath.h"
#include "cdsExcept.h"

#include <iostream>
using namespace std;

class CDSNewtonRaphson {
public:
    int    maxMin;
    int    maxIters;
    bool   verbose;
    double tol;
    double zTol;  // minimization will quit if gradient norm is smaller
    ostream& errStream;

    CDSNewtonRaphson(ostream& errStream=cerr)
      : maxMin( 20 ) , maxIters( 20 ), verbose(0), tol(1e-8), 
        zTol( 1e-8 ), errStream(errStream)
    { }

    template<class CalcB,class CalcZ,class VecType>
    void calc(VecType& x, CalcB calcB, CalcZ calcZ) 
    {
        if ( verbose )
            errStream << "CDSNewtonRaphson: start.\n";
        VecType b = calcB(x);
        double norm = sqrt(abs2(b)) / x.size();
        double zTol2 = CDSMath::sq(zTol * x.size());
        double onorm=norm;
        int  iters=0;
        if ( verbose )
            errStream << "CDSNewtonRaphson: iter: " 
                      << iters << "  norm: " << norm << '\n';
        bool finished=(norm < tol);
        while (!finished) {
            VecType ox = x;
            
            VecType z = calcZ(b);
            x += z;
            
            iters++;

            b = calcB(x);
            norm = sqrt(abs2(b)) / x.size();

            if ( norm > onorm && verbose)
                errStream << "CDSNewtonRaphson: newton-Raphson failed."
                          << " Trying gradient search.\n";

            int mincnt = maxMin;
            z *= -1.0;
            while ( norm > onorm ) {
                if ( mincnt < 1 ) 
                    throw Fail("CDSNewtonRaphson: too many minimization steps taken.\n");
                z *= 0.5;
                if ( abs2(z) < zTol2 )
                    throw Fail("CDSNewtonRaphson: minimization: gradient too small.\n");
                x = ox + z;
                b = calcB(x);
                norm = sqrt(abs2(b)) / x.size();
                if ( verbose )
                    errStream << "CDSNewtonRaphson: iter: " 
                              << iters << "  norm: " << norm << '\n';
                mincnt--;
            }

            onorm = norm;

            if ( verbose )
                errStream << "CDSNewtonRaphson: iter: " 
                          << iters << "  norm: " << norm << '\n';
            
            if (norm < tol)
                finished=1;
            if (iters > maxIters) 
                throw Fail("CDSNewtonRaphson: maxIters exceeded");
        }
    }

    //exception
    struct Fail : public CDS_NAMESPACE(exception) {
        Fail(const char* m) : CDS_NAMESPACE(exception)(m) {}
    };
};

#endif // CDS_NEWTON_RAPHSON_H_
