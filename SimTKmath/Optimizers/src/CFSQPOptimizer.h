#ifndef SimTK_CFSQP_OPTIMIZER_H_
#define SimTK_CFSQP_OPTIMIZER_H_

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
 * Contributors: Eran Guendelman, Frank C. Anderson                           *
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

//=============================================================================
//=============================================================================
/**
 * This class provides methods for finding the optimal controls of a redundant
 * system by applying sequential quadratic programming techniques.  The core
 * algorithm is called fsqp, "Fast Sequential Quadratic Programming".
 */
namespace SimTK { 

class CFSQPOptimizer: public Optimizer::OptimizerRep {
//=============================================================================
// DATA
//=============================================================================
private:
    // PARAMETETERS
    /** Mode. */
    int _mode;
    /** Variable that contains information about the status of an
    optimization.*/
    int _inform;
    /** Value used for infinity or a very large number. */
    double _infinity;
    /** Convergence criterion for nonlinear equality constraints. */
    double _epseqn;
    /** I don't know. */
    double _udelta;

    // SEQUENTIAL INFO
    int _ncsrl;
    int _ncsrn;
    int _nfsr;
    int *_mesh;

    // ALLOCATIONS
    /** Array of performance criteria. */
    double *_p;
    /** Array of constraints. */
    double *_c;
    /** Lagrange multipliers for constraints. */
    double *_lambda;

//=============================================================================
// METHODS
//=============================================================================
public:
    //--------------------------------------------------------------------------
    // CONSTRUCTION
    //--------------------------------------------------------------------------
    virtual ~CFSQPOptimizer();
    CFSQPOptimizer(const OptimizerSystem& sys); 

    static bool isAvailable();

    Real optimize(Vector &results);
    OptimizerRep* clone() const;

    OptimizerAlgorithm getAlgorithm() const
    {   return CFSQP; }

private:
    static void bindToCFSQPLibrary();

    //--------------------------------------------------------------------------
    // STATIC FUNCTIONS USED AS INPUT TO cfsqp()
    //--------------------------------------------------------------------------
    static void
        pFunc(int nparam,int j,double *x,double *pj,void *cd);
    static void
        cFunc(int nparam,int j,double *x,double *cj,void *cd);
    static void
        dpdxFunc(int nparam,int j,double *x,double *dpdx,
        void (*dummy)(int,int,double *,double *,void *),void *cd);
    static void
        dcdxFunc(int nparam,int j,double *x,double *dcdx,
        void (*dummy)(int,int,double *,double *,void *),void *cd);

    //--------------------------------------------------------------------------
    // PRINT
    //--------------------------------------------------------------------------
    static void
        PrintInform(int aInform,std::ostream &aOStream);

    //--------------------------------------------------------------------------
    // CACHING
    //--------------------------------------------------------------------------
    mutable SimTK::Vector _cachedConstraintJacobianParameters;
    mutable SimTK::Matrix _cachedConstraintJacobian;
    mutable SimTK::Vector _cachedConstraintParameters;
    mutable SimTK::Vector _cachedConstraint;

    // Caching support for FSQP (since it likes to query constraints one at a time)
    void clearCache();
public:
    int computeConstraint(const SimTK::Vector &x, const bool new_coefficients, double &c, int ic) const;
    int computeConstraintGradient(const SimTK::Vector &x, const bool new_coefficients, SimTK::Vector &dcdx, int ic) const;

//=============================================================================
};

} // namespace SimTK
//=============================================================================
//=============================================================================

#endif // SimTK_CFSQP_OPTIMIZER_H_
