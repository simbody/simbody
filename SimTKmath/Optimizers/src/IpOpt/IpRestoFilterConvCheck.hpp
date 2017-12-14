// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpRestoFilterConvCheck.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPRESTOFILTERCONVCHECK_HPP__
#define __IPRESTOFILTERCONVCHECK_HPP__

#include "IpOptErrorConvCheck.hpp"
#include "IpFilterLSAcceptor.hpp"

namespace SimTKIpopt
{

  /** Convergence check for the restoration phase as called by the
   *  filter.  This inherits from the OptimalityErrorConvergenceCheck
   *  so that the method for the regular optimality error convergence
   *  criterion can be checked as well.  In addition, this convergence
   *  check returns the CONVERGED message, if the current iteration is
   *  acceptable to the original filter.
   *
   *  Since this object needs to know about the original NLP, it also
   *  inherits from RestoProblemCoupler, so that the restoration phase
   *  object can call the SetObjs method to set the corresponding
   *  pointers before the Initilize for the restoration phase
   *  algorithm is called.
   */
  class RestoFilterConvergenceCheck :
        public OptimalityErrorConvergenceCheck
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    RestoFilterConvergenceCheck();

    /** Default destructor */
    virtual ~RestoFilterConvergenceCheck();
    //@}

    /** Set the object for the original filter line search. Here,
     *  filter_line_search must be the same strategy object to which
     *  the restoration phase object with this object is given.  This
     *  method must be called to finish the definition of the
     *  algorithm, before Initialize is called. */
    void SetOrigFilterLSAcceptor(const FilterLSAcceptor& orig_filter_ls_acceptor);

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** overloaded from ConvergenceCheck */
    virtual ConvergenceStatus CheckConvergence(bool call_intermediate_callback = true) override;

    /** Methods used by IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}
  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    RestoFilterConvergenceCheck(const RestoFilterConvergenceCheck&);

    /** Overloaded Equals Operator */
    void operator=(const RestoFilterConvergenceCheck&);
    //@}

    /** @name Algorithmic parameters */
    //@{
    /** Fraction of required reduction in infeasibility before problem
     *  is considered to be solved. */
    Number kappa_resto_;
    /** Maximum number of iterations in restoration phase */
    Index maximum_iters_;
    //@}

    /** Flag indicating that this is the first call.  We don't want to
     *  leave the restoration phase without taking at least one step,
     *  so this flag is used to ensure this. */
    bool first_resto_iter_;

    /** Strategy object for the filter line search method for the
     *  original NLP.  CAREFUL: We must not hold on to this object
     *  with a SmartPtr, because have otherwise circular references
     *  that prevent the destructor of the line search object to be
     *  called! */
    const FilterLSAcceptor* orig_filter_ls_acceptor_;
  };

} // namespace Ipopt

#endif
