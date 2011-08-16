#ifndef _SimTK_LBFGSB_OPTIMIZER_H_
#define _SimTK_LBFGSB_OPTIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-10 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
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

#include "simmath/internal/common.h"
#include "OptimizerRep.h"

namespace SimTK {

class LBFGSBOptimizer: public Optimizer::OptimizerRep {
public:
    ~LBFGSBOptimizer() {
        delete [] nbd;
    }

    LBFGSBOptimizer(const OptimizerSystem& sys); 

    Real optimize(  Vector &results );
    OptimizerRep* clone() const;

    int setulb_(int *n, int *m, Real *x, Real *l,
        Real *u, int *nbd, Real *f, Real *g,
        Real *factr, Real *pgtol, Real *wa, int *iwa,
        char *task, int *iprint, char *csave, bool *lsave,
        int *isave, Real *dsave, long task_len, long csave_len);

private:
    double      factr;
    int         iprint[3];
    int         *nbd;
};

} // namespace SimTK
#endif //_SimTK_LBFGSB_OPTIMIZER_H_

