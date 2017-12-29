// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIterationOutput.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter, Carl Laird       IBM    2004-09-27

#ifndef __IPITERATIONOUTPUT_HPP__
#define __IPITERATIONOUTPUT_HPP__

#include "IpAlgStrategy.hpp"
#include "IpIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace SimTKIpopt
{

  /** Base class for objects that do the output summary per iteration.
   */
  class IterationOutput: public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    IterationOutput()
    {}

    /** Default destructor */
    virtual ~IterationOutput()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Method to do all the summary output per iteration.  This
     *  include the one-line summary output as well as writing the
     *  details about the iterates if desired */
    virtual void WriteOutput() = 0;

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    IterationOutput(const IterationOutput&);

    /** Overloaded Equals Operator */
    void operator=(const IterationOutput&);
    //@}

  };

} // namespace Ipopt

#endif
