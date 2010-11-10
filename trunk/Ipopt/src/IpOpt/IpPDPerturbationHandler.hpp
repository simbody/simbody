// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpPDPerturbationHandler.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2005-08-04

#ifndef __IPPDPERTURBATIONHANDLER_HPP__
#define __IPPDPERTURBATIONHANDLER_HPP__

#include "IpAlgStrategy.hpp"

namespace Ipopt
{

  /** Class for handling the perturbation factors delta_x, delta_s,
   *  delta_c, and delta_d in the primal dual system.  This class is
   *  used by the PDFullSpaceSolver to handle the cases where the
   *  primal-dual system is singular or has the wrong inertia.  The
   *  perturbation factors are obtained based on simple heuristics,
   *  taking into account the size of previous perturbations.
   */
  class PDPerturbationHandler: public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    PDPerturbationHandler();
    /** Default destructor */
    virtual ~PDPerturbationHandler()
    {}
    //@}

    /* overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** This method must be called for each new matrix, and before any
     *  other method for generating perturbation factors.  Usually,
     *  the returned perturbation factors are zero, but if the system
     *  is thought to be structurally singular, they might be
     *  positive.  If the return value is false, no suitable
     *  perturbation could be found. */
    bool ConsiderNewSystem(Number& delta_x, Number& delta_s,
                           Number& delta_c, Number& delta_d);

    /** This method returns pertubation factors for the case when the
     *  most recent factorization resulted in a singular matrix. If
     *  the return value is false, no suitable perturbation could be
     *  found. */
    bool PerturbForSingularity(Number& delta_x, Number& delta_s,
                               Number& delta_c, Number& delta_d);

    /** This method returns pertubation factors for the case when the
     *  most recent factorization resulted in a matrix with an
     *  incorrect number of negative eigenvalues. If the return value
     *  is false, no suitable perturbation could be found. */
    bool PerturbForWrongInertia(Number& delta_x, Number& delta_s,
                                Number& delta_c, Number& delta_d);

    /** Just return the perturbation values that have been determined
     *  most recently */
    void CurrentPerturbation(Number& delta_x, Number& delta_s,
                             Number& delta_c, Number& delta_d);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    PDPerturbationHandler(const PDPerturbationHandler&);

    /** Overloaded Equals Operator */
    void operator=(const PDPerturbationHandler&);
    //@}

    /** @name Size of the most recent non-zero perturbation. */
    //@{
    /** The last nonzero value for delta_x */
    Number delta_x_last_;
    /** The last nonzero value for delta_s */
    Number delta_s_last_;
    /** The last nonzero value for delta_c */
    Number delta_c_last_;
    /** The last nonzero value for delta_d */
    Number delta_d_last_;
    //@}

    /** @name Size of the most recently suggested perturbation for the
     *  current matrix. */
    //@{
    /** The current value for delta_x */
    Number delta_x_curr_;
    /** The current value for delta_s */
    Number delta_s_curr_;
    /** The current value for delta_c */
    Number delta_c_curr_;
    /** The current value for delta_d */
    Number delta_d_curr_;
    //@}

    /** Flag indicating if for the given matrix the perturb for wrong
     *  inertia method has already been called. */
    bool get_deltas_for_wrong_inertia_called_;

    /** @name Handling structural degeneracy */
    //@{
    /** Type for degeneracy flags */
    enum DegenType {
      NOT_YET_DETERMINED,
      NOT_DEGENERATE,
      DEGENERATE
    };

    /** Flag indicating whether the reduced Hessian matrix is thought
     *  to be structurally singular. */
    DegenType hess_degenerate_;

    /** Flag indicating whether the Jacobian of the constraints is
     *  thought to be structurally rank-deficient. */
    DegenType jac_degenerate_;

    /** Flag counting matrices in which degeneracy was observed in the
     *  first successive iterations.  -1 means that there was a
     *  non-degenerate (unperturbed) matrix at some point. */
    Index degen_iters_;

    /** Status of current trial configuration */
    enum TrialStatus {
      NO_TEST,
      TEST_DELTA_C_EQ_0_DELTA_X_EQ_0,
      TEST_DELTA_C_GT_0_DELTA_X_EQ_0,
      TEST_DELTA_C_EQ_0_DELTA_X_GT_0,
      TEST_DELTA_C_GT_0_DELTA_X_GT_0
    };

    /** Current status */
    TrialStatus test_status_;
    //@}

    /** @name Algorithmic parameters. */
    //@{
    /** Maximal perturbation for x and s. */
    Number delta_xs_max_;
    /** Smallest possible perturbation for x and s. */
    Number delta_xs_min_;
    /** Increase factor for delta_xs for first required perturbation. */
    Number delta_xs_first_inc_fact_;
    /** Increase factor for delta_xs for later perturbations. */
    Number delta_xs_inc_fact_;
    /** Decrease factor for delta_xs for later perturbations. */
    Number delta_xs_dec_fact_;
    /** Very first trial value for delta_xs perturbation. */
    Number delta_xs_init_;
    /** Size of perturbation for c and d blocks. */
    Number delta_cd_val_;
    /** Exponent on mu in formula for of perturbation for c and d blocks. */
    Number delta_cd_exp_;
    /** Flag indicating whether the new values are based on the
     *  perturbations in the last iteration or in the more recent
     *  iteration in which a perturbation was done. */
    bool reset_last_;
    /** Required number of iterations for degeneracy conclusions. */
    Index degen_iters_max_;
    //@}

    /** @name Auxilliary methods */
    //@{
    /** Internal version of PerturbForWrongInertia with the
     *  difference, that finalize_test is not called.  Returns false
     *  if the delta_x and delta_s parameters become too large. */
    bool get_deltas_for_wrong_inertia(Number& delta_x, Number& delta_s,
                                      Number& delta_c, Number& delta_d);

    /** This method is call whenever a matrix had been factorization
     *  and is not singular.  In here, we can evaluate the outcome of
     *  the deneracy test heuristics. */
    void finalize_test();
    /** Compute perturbation value for constraints */
    Number delta_cd();
    //@}

  };

} // namespace Ipopt

#endif
