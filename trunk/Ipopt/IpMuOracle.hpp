// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMuOracle.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPMUORACLE_HPP__
#define __IPMUORACLE_HPP__

#include "IpAlgStrategy.hpp"

namespace Ipopt
{

  /** Abstract Base Class for classes that are able to compute a
   *  suggested value of the barrier parameter that can be used
   *  as an oracle in the NonmontoneMuUpdate class.
   */
  class MuOracle : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    MuOracle()
    {}
    ;
    /** Default destructor */
    virtual ~MuOracle()
    {}
    ;
    //@}

    /** Initialize method - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Method for computing the value of the barrier parameter that
     *  could be used in the current iteration.  Here, mu_min and
     *  mu_max are the lower and upper bounds on acceptable values for
     *  the barrier parameter.  The new value of mu is returned in
     *  new_mu, and the method returns false if a new value could not
     *  be determined (e.g., because the linear system could not be
     *  solved for a predictor step).
     */
    virtual bool CalculateMu(Number mu_min, Number mu_max, Number& new_mu) = 0;

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
    MuOracle(const MuOracle&);

    /** Overloaded Equals Operator */
    void operator=(const MuOracle&);
    //@}

  };

} // namespace Ipopt

#endif
