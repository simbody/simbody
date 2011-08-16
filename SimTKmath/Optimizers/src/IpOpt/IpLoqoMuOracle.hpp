// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLoqoMuOracle.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPLOQOMUORACLE_HPP__
#define __IPLOQOMUORACLE_HPP__

#include "IpMuOracle.hpp"

namespace Ipopt
{

  /** Implementation of the LOQO formula for computing the
   *  barrier parameter.
   */
  class LoqoMuOracle : public MuOracle
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    LoqoMuOracle();
    /** Default destructor */
    virtual ~LoqoMuOracle();
    //@}

    /** Initialize method - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the value of the barrier parameter that
     *  could be used in the current iteration (using the LOQO formula).
     */
    virtual bool CalculateMu(Number mu_min, Number mu_max, Number& new_mu);

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
    LoqoMuOracle(const LoqoMuOracle&);

    /** Overloaded Equals Operator */
    void operator=(const LoqoMuOracle&);
    //@}

  };

} // namespace Ipopt

#endif
