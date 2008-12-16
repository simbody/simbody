#ifndef _SimTK_CFSQP_OPTIMIZER_H_
#define _SimTK_CFSQP_OPTIMIZER_H_

/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
 * Contributors: Eran Guendelman, Frank C. Anderson
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "SimTKcommon.h"

#include "simmath/internal/common.h"

#include "OptimizerRep.h"

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

#endif // _SimTK_CFSQP_OPTIMIZER_H_
