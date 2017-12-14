// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptApplication.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTAPPLICATION_HPP__
#define __IPIPOPTAPPLICATION_HPP__

#include <iostream>

#include "IpJournalist.hpp"
#include "IpTNLP.hpp"
#include "IpNLP.hpp"
/* Return codes for the Optimize call for an application */
#include "IpReturnCodes.hpp"

namespace SimTKIpopt
{
  DECLARE_STD_EXCEPTION(IPOPT_APPLICATION_ERROR);

  /* forward declarations */
  class IpoptAlgorithm;
  class IpoptNLP;
  class IpoptData;
  class IpoptCalculatedQuantities;
  class AlgorithmBuilder;
  class RegisteredOptions;
  class OptionsList;
  class SolveStatistics;

  /** This is the main application class for making calls to Ipopt. */
  class IpoptApplication : public ReferencedObject
  {
  public:
    IpoptApplication(bool create_console_out = true);

    virtual ~IpoptApplication();

    /** Initialize method. This method reads the params file and initializes
     *  the journalists. You should call this method at some point before the 
     *  first optimize call. Note: you can skip the processing of a params
     *  file by setting params_file to ""
     */
    void Initialize(std::string params_file = "ipopt.opt");
    void Initialize(std::istream& is);

    /**@name Solve methods */
    //@{
    /** Solve a problem that inherits from TNLP */
    ApplicationReturnStatus OptimizeTNLP(const SmartPtr<TNLP>& tnlp);

    /** Solve a problem that inherits from NLP */
    ApplicationReturnStatus OptimizeNLP(const SmartPtr<NLP>& nlp, SmartPtr<AlgorithmBuilder> alg_builder=NULL);

    /** Solve a problem (that inherits from TNLP) for a repeated time.
     *  The OptimizeTNLP method must have been called before.  The
     *  TNLP must be the same object, and the structure (number of
     *  variables and constraints and position of nonzeros in Jacobian
     *  and Hessian must be the same). */
    ApplicationReturnStatus ReOptimizeTNLP(const SmartPtr<TNLP>& tnlp);

    /** Solve a problem (that inherits from NLP) for a repeated time.
     *  The OptimizeNLP method must have been called before.  The
     *  NLP must be the same object, and the structure (number of
     *  variables and constraints and position of nonzeros in Jacobian
     *  and Hessian must be the same). */
    ApplicationReturnStatus ReOptimizeNLP(const SmartPtr<NLP>& nlp);
    //@}

    /** Method for opening an output file with given print_level.
     *  Returns false if there was a problem. */
    bool OpenOutputFile(std::string file_name, EJournalLevel print_level);

    /**@name Accessor methods */
    //@{
    /** Get the Journalist for printing output */
    SmartPtr<Journalist> Jnlst()
    {
      return jnlst_;
    }

    /** Get a pointer to RegisteredOptions object to
     *  add new options */
    SmartPtr<RegisteredOptions> RegOptions()
    {
      return reg_options_;
    }

    /** Get the options list for setting options */
    SmartPtr<OptionsList> Options()
    {
      return options_;
    }

    /** Get the options list for setting options (const version) */
    SmartPtr<const OptionsList> Options() const
    {
      return ConstPtr(options_);
    }

    /** Get the object with the statistics about the most recent
     *  optimization run. */
    SmartPtr<SolveStatistics> Statistics();
    //@}

    /** @name Methods for IpoptTypeInfo */
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
    /** Default Constructor */
    // IpoptApplication();

    /** Copy Constructor */
    IpoptApplication(const IpoptApplication&);

    /** Overloaded Equals Operator */
    void operator=(const IpoptApplication&);
    //@}

    /** Method to register all the options */
    void RegisterAllOptions(const SmartPtr<RegisteredOptions>& roptions);

    /** Method for the actual optimize call of the Ipopt algorithm.
     *  This is used both for Optimize and ReOptimize */
    ApplicationReturnStatus call_optimize();

    /**@name Variables that customize the application behavior */
    //@{
    /** Decide whether or not the ipopt.opt file should be read */
    bool read_params_dat_;
    //@}

    /** Journalist for reporting output */
    SmartPtr<Journalist> jnlst_;

    /** RegisteredOptions */
    SmartPtr<RegisteredOptions> reg_options_;

    /** OptionsList used for the application */
    SmartPtr<OptionsList> options_;

    /** Object for storing statistics about the most recent
     *  optimization run. */
    SmartPtr<SolveStatistics> statistics_;

    /** Object with the algorithm sceleton.
     */
    SmartPtr<IpoptAlgorithm> alg_;

    /** IpoptNLP Object for the NLP.  We keep this around for a
     *  ReOptimize warm start. */
    SmartPtr<IpoptNLP> ip_nlp_;

    /** IpoptData Object for the NLP.  We keep this around for a
     *  ReOptimize warm start.
     */
    SmartPtr<IpoptData> ip_data_;

    /** IpoptCalculatedQuantities Object for the NLP.  We keep this
     *  around for a ReOptimize warm start.
     */
    SmartPtr<IpoptCalculatedQuantities> ip_cq_;

    /** Pointer to the TNLPAdapter used to convert the TNLP to an NLP.
     *  We keep this around for the ReOptimizerTNLP call. */
    SmartPtr<NLP> nlp_adapter_;
  };

} // namespace Ipopt

#endif
