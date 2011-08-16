// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpRestoIterateInitializer.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-09-24

#ifndef __IPRESTOITERATEINITIALIZER_HPP__
#define __IPRESTOITERATEINITIALIZER_HPP__

#include "IpIterateInitializer.hpp"
#include "IpEqMultCalculator.hpp"

namespace Ipopt
{
  /** Class implementing the default initialization procedure (based
   *  on user options) for the iterates.  It is used at the very
   *  beginning of the optimization for determine the starting point
   *  for all variables.
   */
  class RestoIterateInitializer: public IterateInitializer
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  If eq_mult_calculator is not NULL, it will be
     *  used to compute the initial values for equality constraint
     *  multipliers. */
    RestoIterateInitializer
    (const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator);

    /** Default destructor */
    virtual ~RestoIterateInitializer()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Compute the initial iterates and set the into the curr field
     *  of the ip_data object. */
    virtual bool SetInitialIterates();

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
    RestoIterateInitializer();

    /** Copy Constructor */
    RestoIterateInitializer(const RestoIterateInitializer&);

    /** Overloaded Equals Operator */
    void operator=(const RestoIterateInitializer&);
    //@}

    /**@name Parameters for bumping x0 */
    //@{
    /** If max-norm of the initial equality constraint multiplier
     *  estimate is larger than this, the initial y_* variables are
     *  set to zero. */
    Number constr_mult_init_max_;
    //@}

    /** object to be used for the initialization of the equality
     *  constraint multipliers. */
    SmartPtr<EqMultiplierCalculator> resto_eq_mult_calculator_;

    /** @name Auxilliary functions */
    //@{
    /** Method for solving the quadratic vector equation v^2 + 2a*v -
    b = 0 */
    void solve_quadratic(const Vector& a, const Vector& b, Vector& v);
    //@}
  };

} // namespace Ipopt

#endif
