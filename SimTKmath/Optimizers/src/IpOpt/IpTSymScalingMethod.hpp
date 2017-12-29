// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpTSymScalingMethod.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPTSYMSCALINGMETHOD_HPP__
#define __IPTSYMSCALINGMETHOD_HPP__

#include "IpUtils.hpp"
#include "IpAlgStrategy.hpp"

namespace SimTKIpopt
{

  DECLARE_STD_EXCEPTION(ERROR_IN_LINEAR_SCALING_METHOD);

  /** Base class for the method for computing scaling factors for symmetric
   *  matrices in triplet format.
   */
  class TSymScalingMethod: public AlgorithmStrategyObject
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    TSymScalingMethod()
    {}

    ~TSymScalingMethod()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Method for computing the symmetric scaling factors, given the
     *  symmtric matrix in triplet (MA27) format. */
    virtual bool ComputeSymTScalingFactors(Index n,
                                           Index nnz,
                                           const Index* airn,
                                           const Index* ajcn,
                                           const Number* a,
                                           Number* scaling_factors) = 0;
  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    TSymScalingMethod(const TSymScalingMethod&);

    /** Overloaded Equals Operator */
    void operator=(const TSymScalingMethod&);
  };

} // namespace Ipopt

#endif
