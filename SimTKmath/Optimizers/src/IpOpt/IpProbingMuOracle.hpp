// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpProbingMuOracle.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPPROBINGMUORACLE_HPP__
#define __IPPROBINGMUORACLE_HPP__

#include "IpMuOracle.hpp"
#include "IpPDSystemSolver.hpp"

namespace SimTKIpopt
{

  /** Implementation of the probing strategy for computing the
   *  barrier parameter.
   */
  class ProbingMuOracle : public MuOracle
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    ProbingMuOracle(const SmartPtr<PDSystemSolver>& pd_solver);
    /** Default destructor */
    virtual ~ProbingMuOracle();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** Method for computing the value of the barrier parameter that
     *  could be used in the current iteration (using Mehrotra's
     *  probing heuristic).
     */
    virtual bool CalculateMu(Number mu_min, Number mu_max, Number& new_mu) override;

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

    /** Default Constructor */
    ProbingMuOracle();
    /** Copy Constructor */
    ProbingMuOracle(const ProbingMuOracle&);

    /** Overloaded Equals Operator */
    void operator=(const ProbingMuOracle&);
    //@}

    /** Pointer to the object that should be used to solve the
     *  primal-dual system.
     */
    SmartPtr<PDSystemSolver> pd_solver_;

    /** Auxilliary function for computing the average complementarity
     *  at a point, given step sizes and step
     */
    Number CalculateAffineMu(Number alpha_primal,
                             Number alpha_dual,
                             const IteratesVector& step);

    /** @name Algorithmic parameters */
    //@{
    /** safeguarding upper bound on centering parameter sigma */
    Number sigma_max_;
    //@}
  };

} // namespace Ipopt

#endif
