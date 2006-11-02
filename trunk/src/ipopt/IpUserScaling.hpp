// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpUserScaling.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-25

#ifndef __IPUSERSCALING_HPP__
#define __IPUSERSCALING_HPP__

#include "IpNLPScaling.hpp"
#include "IpNLP.hpp"

namespace Ipopt
{
  /** This class does problem scaling by getting scaling parameters
   *  from the user (through the NLP interface).
   */
  class UserScaling : public StandardScalingBase
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    UserScaling(const SmartPtr<const NLP>& nlp)
        :
        StandardScalingBase(),
        nlp_(nlp)
    {}

    /** Default destructor */
    virtual ~UserScaling()
    {}
    //@}

  protected:
    virtual void DetermineScalingParametersImpl(
      const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      const SmartPtr<const MatrixSpace> jac_c_space,
      const SmartPtr<const MatrixSpace> jac_d_space,
      const SmartPtr<const SymMatrixSpace> h_space,
      Number& df,
      SmartPtr<Vector>& dx,
      SmartPtr<Vector>& dc,
      SmartPtr<Vector>& dd);

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
    UserScaling(const UserScaling&);

    /** Overloaded Equals Operator */
    void operator=(const UserScaling&);
    //@}

    /** pointer to the NLP to get scaling parameters */
    SmartPtr<const NLP> nlp_;
  };
} // namespace Ipopt
#endif
