// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpWarmStartIterateInitializer.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2005-04-01

#ifndef __IPWARMSTARTITERATEINITIALIZER_HPP__
#define __IPWARMSTARTITERATEINITIALIZER_HPP__

#include "IpIterateInitializer.hpp"
#include "IpEqMultCalculator.hpp"

namespace SimTKIpopt
{

  /** Class implementing an initialization procedure for warm starts.
   */
  class WarmStartIterateInitializer: public IterateInitializer
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor. */
    WarmStartIterateInitializer();

    /** Default destructor */
    virtual ~WarmStartIterateInitializer()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** Compute the initial iterates and set the into the curr field
     *  of the ip_data object. */
    virtual bool SetInitialIterates() override;

    /** Methods used by IpoptType */
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
    WarmStartIterateInitializer(const WarmStartIterateInitializer&);

    /** Overloaded Equals Operator */
    void operator=(const WarmStartIterateInitializer&);
    //@}

    /**@name Algorithmic Parameters */
    //@{
    /** Parameters for bumping x0 in warm start mode */
    Number warm_start_bound_push_;
    /** Parameters for bumping x0 in warm start mode */
    Number warm_start_bound_frac_;
    /** Parameters for bumping initial bound multipliers */
    Number warm_start_mult_bound_push_;
    /** Maximal size of entries in bound and equality constraint
     *  multipliers in magnitute.  If chosen less of equal to zero, no
     *  upper limit is imposed.  Otherwise, the entries exceeding the
     *  given limit are set to the value closest to the limit. */
    Number warm_start_mult_init_max_;
    /** Target values for the barrier parameter in warm start option.
     */
    Number warm_start_target_mu_;
    /** Indicator for which method in the NLP should be used to get
     *  the warm start  */
    bool warm_start_entire_iterate_;
    //@}

    /** @name Auxilliary functions */
    //@{
    void process_target_mu(Number factor,
                           const Vector& curr_vars,
                           const Vector& curr_slacks,
                           const Vector& curr_mults,
                           const Matrix& P,
                           SmartPtr<const Vector>& ret_vars,
                           SmartPtr<const Vector>& ret_mults);

    void adapt_to_target_mu(Vector& new_s,
                            Vector& new_z,
                            Number target_mu);
    //@}
  };

} // namespace Ipopt

#endif
