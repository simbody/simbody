// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpRestoPhase.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPRESTOPHASE_HPP__
#define __IPRESTOPHASE_HPP__

#include "IpAlgStrategy.hpp"
#include "IpIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace SimTKIpopt
{

  /** @name Exceptions */
  //@{
  /** Exception RESTORATION_FAILED for all trouble with the
   *  restoration phase.
   */
  DECLARE_STD_EXCEPTION(RESTORATION_CONVERGED_TO_FEASIBLE_POINT);
  DECLARE_STD_EXCEPTION(RESTORATION_FAILED);
  DECLARE_STD_EXCEPTION(RESTORATION_MAXITER_EXCEEDED);
  DECLARE_STD_EXCEPTION(RESTORATION_USER_STOP);
  //@}

  /** Base class for different restoration phases.  The restoration
   *  phase is part of the FilterLineSearch. */
  class RestorationPhase : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    RestorationPhase()
    {}
    /** Default Destructor */
    virtual ~RestorationPhase()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Method called to perform restoration for the filter line
     *  search method. */
    virtual bool PerformRestoration() = 0;

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    RestorationPhase(const RestorationPhase&);

    /** Overloaded Equals Operator */
    void operator=(const RestorationPhase&);
    //@}
  };

} // namespace Ipopt

#endif
