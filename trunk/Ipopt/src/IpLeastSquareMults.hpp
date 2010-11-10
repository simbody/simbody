// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLeastSquareMults.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-09-23

#ifndef __IPLEASTSQUAREMULTS_HPP__
#define __IPLEASTSQUAREMULTS_HPP__

#include "IpAugSystemSolver.hpp"
#include "IpEqMultCalculator.hpp"

namespace Ipopt
{

  /** Class for calculator for the least-square equality constraint
   *  multipliers.  The Calculate method of this class computes the
   *  least-square estimate for the y_c and y_d multiplers, based on
   *  the current values of the gradient of the Lagrangian.
   */
  class LeastSquareMultipliers: public EqMultiplierCalculator
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  It needs to be given the strategy object for
     *  solving the augmented system. */
    LeastSquareMultipliers(AugSystemSolver& augSysSolver);
    /** Default destructor */
    virtual ~LeastSquareMultipliers()
    {}
    //@}

    /* overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** This method computes the least-square estimates for y_c and
     *  y_d at the current point.  The return value is false, if the
     *  least square system could not be solved (the linear system is
     *  singular). */
    virtual bool CalculateMultipliers(Vector& y_c,
                                      Vector& y_d);

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
    LeastSquareMultipliers();

    /** Copy Constructor */
    LeastSquareMultipliers(const LeastSquareMultipliers&);

    /** Overloaded Equals Operator */
    void operator=(const LeastSquareMultipliers&);
    //@}

    /** Pointer for the augmented system solver to be used for solving
     *  the linear system */
    SmartPtr<AugSystemSolver> augsyssolver_;
  };

} // namespace Ipopt

#endif
