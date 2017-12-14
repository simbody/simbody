// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpBacktrackingLSAcceptor.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.hpp

#ifndef __IPBACKTRACKINGLSACCEPTOR_HPP__
#define __IPBACKTRACKINGLSACCEPTOR_HPP__

#include "IpAlgStrategy.hpp"

namespace SimTKIpopt
{

  /** Base class for backtracking line search acceptors.
   */
  class BacktrackingLSAcceptor : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor. */
    BacktrackingLSAcceptor()
    {}

    /** Default destructor */
    virtual ~BacktrackingLSAcceptor()
    {}
    //@}

    /** InitializeImpl - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix)=0;

    /** Reset the acceptor.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be called if
     *  the barrier parameter is changed.
     */
    virtual void Reset()=0;

    /** Initialization for the next line search.  The flag in_watchdog
     *  indicates if we are currently in an active watchdog
     *  procedure. */
    virtual void InitThisLineSearch(bool in_watchdog)=0;

    /** Method that is called before the restoration phase is called.
     *  Here, we can set up things that are required in the
     *  termination test for the restoration phase, such as augmenting
     *  a filter. */
    virtual void PrepareRestoPhaseStart()=0;

    /** Method returning the lower bound on the trial step sizes.  If
     *  the backtracking procedure encounters a trial step size below
     *  this value after the first trial set, it swtiches to the
     *  (soft) restoration phase. */
    virtual Number CalculateAlphaMin()=0;

    /** Method for checking if current trial point is acceptable.  It
     *  is assumed that the delta information in ip_data is the search
     *  direction used in criteria.  The primal trial point has to be
     *  set before the call.  alpha_primal is the step size which is
     *  to be used for the test; if it is zero, then this method is
     *  called during the soft restoration phase.
     */
    virtual bool CheckAcceptabilityOfTrialPoint(Number alpha_primal)=0;

    /** Try a second order correction for the constraints.  If the
     *  first trial step (with incoming alpha_primal) has been reject,
     *  this tries second order corrections, e.g., for the
     *  constraints.  Here, alpha_primal_test is the step size that
     *  has to be used in the filter acceptance tests.  On output
     *  actual_delta_ has been set to the step including the
     *  second order correction if it has been accepted, otherwise it
     *  is unchanged.  If the SOC step has been accepted, alpha_primal
     *  has the fraction-to-the-boundary value for the SOC step on output.
     *  The return value is true, if a SOC step has been accepted.
     */
    virtual bool TrySecondOrderCorrection(Number alpha_primal_test,
                                          Number& alpha_primal,
                                          SmartPtr<IteratesVector>& actual_delta)=0;

    /** Try higher order corrector (for fast local convergence).  In
     *  contrast to a second order correction step, which tries to
     *  make an unacceptable point acceptable by improving constraint
     *  violation, this corrector step is tried even if the regular
     *  primal-dual step is acceptable.
     */
    virtual bool TryCorrector(Number alpha_primal_test,
                              Number& alpha_primal,
                              SmartPtr<IteratesVector>& actual_delta)=0;

    /** Method for ending the current line search.  When it is called,
     *  the internal data should be updates, e.g., the filter might be
     *  augmented.  alpha_primal_test is the value of alpha that has
     *  been used for in the acceptance test earlier.  Return value is
     *  a character for the info_alpha_primal_char field in IpData. */
    virtual char UpdateForNextIteration(Number alpha_primal_test)=0;

    /** Method for setting internal data if the watchdog procedure is
     *  started. */
    virtual void StartWatchDog()=0;

    /** Method for setting internal data if the watchdog procedure is
     *  stopped. */
    virtual void StopWatchDog()=0;

    /** Methods for OptionsList */
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
    BacktrackingLSAcceptor(const BacktrackingLSAcceptor&);

    /** Overloaded Equals Operator */
    void operator=(const BacktrackingLSAcceptor&);
    //@}
  };

} // namespace Ipopt

#endif
