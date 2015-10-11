#ifndef SimTK_SIMMATH_LBFGS_OPTIMIZER_H_
#define SimTK_SIMMATH_LBFGS_OPTIMIZER_H_

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

#include "SimTKcommon.h"

#include "simmath/internal/common.h"

#include "simmath/internal/OptimizerRep.h"

#include <iostream>

namespace SimTK {


class LBFGSOptimizer: public Optimizer::OptimizerRep {
public:
    ~LBFGSOptimizer() { }

    LBFGSOptimizer(const OptimizerSystem& sys);

    Real optimize(  SimTK::Vector &results ) override;
    OptimizerRep* clone() const override;

    OptimizerAlgorithm getAlgorithm() const override
    {   return LBFGS; }

    private:
    int         iprint[3];
    Real        xtol;
    void lbfgs_( int n, int m, Real *x, Real *f, int *iprint,  Real *eps, Real *xtol );
    void mcsrch_(int *n, Real *x, Real *f, Real *g, Real *s, Real *stp,
                 Real *ftol, Real *xtol, int *maxfev, int *info, int *nfev, Real *wa);

    void setXtol( double );
};

} // namespace SimTK

#endif // SimTK_SIMMATH_LBFGS_OPTIMIZER_H_

