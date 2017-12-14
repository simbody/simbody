// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOrigIterationOutput.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter, Carl Laird       IBM    2004-09-27

#ifndef __IPORIGITERATIONOUTPUT_HPP__
#define __IPORIGITERATIONOUTPUT_HPP__

#include "IpIterationOutput.hpp"

namespace SimTKIpopt
{

  /** Class for the iteration summary output for the original NLP.
   */
  class OrigIterationOutput: public IterationOutput
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    OrigIterationOutput();

    /** Default destructor */
    virtual ~OrigIterationOutput();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** Method to do all the summary output per iteration.  This
     *  include the one-line summary output as well as writing the
     *  details about the iterates if desired */
    virtual void WriteOutput() override;

    /** Methods for OptionsList */
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
    OrigIterationOutput(const OrigIterationOutput&);

    /** Overloaded Equals Operator */
    void operator=(const OrigIterationOutput&);
    //@}

    /** Flag indicating weather info string should be printed at end
     *  of iteration summary line. */
    bool print_info_string_;
  };

} // namespace Ipopt

#endif
