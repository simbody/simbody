#ifndef _SimTK_LBFGSB_OPTIMIZER_H_
#define _SimTK_LBFGSB_OPTIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
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

#include "simmath/internal/common.h"
#include "simmath/internal/OptimizerRep.h"

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
    Real        factr;
    int         iprint[3];
    int         *nbd;
};

} // namespace SimTK
#endif //_SimTK_LBFGSB_OPTIMIZER_H_

