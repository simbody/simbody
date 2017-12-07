// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpRestoRestoPhase.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-02-11

#ifndef __IPRESTORESTOPHASE_HPP__
#define __IPRESTORESTOPHASE_HPP__

#include "IpRestoPhase.hpp"
#include "IpIpoptAlg.hpp"
#include "IpEqMultCalculator.hpp"

namespace SimTKIpopt
{

  /** Recursive Restoration Phase for the.MinC_1NrmRestorationPhase.
   *  This procedure chooses the n and p variables in the
   *  MinC_1NrmRestorationPhase problem formulation by treating the
   *  problem as separable (assuming that the x and s variables are
   *  fixed).
   */
  class RestoRestorationPhase : public RestorationPhase
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor. */
    RestoRestorationPhase();

    /** Default destructor */
    virtual ~RestoRestorationPhase();
    //@}

    /** Overloaded from AlgorithmStrategy case class */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

  protected:
    /** Overloaded method from RestorationPhase. */
    virtual bool PerformRestoration() override;

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    RestoRestorationPhase(const RestoRestorationPhase&);

    /** Overloaded Equals Operator */
    void operator=(const RestoRestorationPhase&);
    //@}

    /** @name Auxilliary methods */
    //@{
    /** Method for solving the quadratic vector equation v^2 + 2a*v -
    b = 0 */
    void solve_quadratic(const Vector& a, const Vector& b, Vector& v);
    //@}
  };

} // namespace Ipopt

#endif
