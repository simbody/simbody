// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpTimingStatistics.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter               IBM    2005-09-19

#ifndef __IPTIMINGSTATISTICS_HPP__
#define __IPTIMINGSTATISTICS_HPP__

#include "IpReferenced.hpp"
#include "IpJournalist.hpp"
#include "IpTimedTask.hpp"

#ifdef HAVE_CTIME
# include <ctime>
#else
# ifdef HAVE_TIME_H
#  include <time.h>
# else
#  error "don't have header file for time"
# endif
#endif

// The following lines are copied from CoinTime.hpp
// We should probably make some more tests here
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#else
// MacOS-X and FreeBSD needs sys/time.h
# if defined(__MACH__) || defined (__FreeBSD__)
#  include <sys/time.h>
# endif
# if !defined(__MSVCRT__)
#  include <sys/resource.h>
# endif
#endif

namespace Ipopt
{
  /** This class collects all timing statistics for Ipopt.
   */
  class TimingStatistics : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor. */
    TimingStatistics()
    {}

    /** Default destructor */
    virtual ~TimingStatistics()
    {}
    //@}

    /** Method for resetting all times. */
    void ResetTimes();

    /** Method for printing all timing information */
    void PrintAllTimingStatistics(Journalist& jnlst,
                                  EJournalLevel level,
                                  EJournalCategory category) const;

    /**@name Accessor methods to all timed tasks. */
    //@{
    TimedTask& OverallAlgorithm()
    {
      return OverallAlgorithm_;
    }
    TimedTask& PrintProblemStatistics()
    {
      return PrintProblemStatistics_;
    }
    TimedTask& InitializeIterates()
    {
      return InitializeIterates_;
    }
    TimedTask& UpdateHessian()
    {
      return UpdateHessian_;
    }
    TimedTask& OutputIteration()
    {
      return OutputIteration_;
    }
    TimedTask& UpdateBarrierParameter()
    {
      return UpdateBarrierParameter_;
    }
    TimedTask& ComputeSearchDirection()
    {
      return ComputeSearchDirection_;
    }
    TimedTask& ComputeAcceptableTrialPoint()
    {
      return ComputeAcceptableTrialPoint_;
    }
    TimedTask& AcceptTrialPoint()
    {
      return AcceptTrialPoint_;
    }
    TimedTask& CheckConvergence()
    {
      return CheckConvergence_;
    }

    TimedTask& PDSystemSolverTotal()
    {
      return PDSystemSolverTotal_;
    }
    TimedTask& PDSystemSolverSolveOnce()
    {
      return PDSystemSolverSolveOnce_;
    }
    TimedTask& ComputeResiduals()
    {
      return ComputeResiduals_;
    }
    TimedTask& LinearSystemScaling()
    {
      return LinearSystemScaling_;
    }
    TimedTask& LinearSystemSymbolicFactorization()
    {
      return LinearSystemSymbolicFactorization_;
    }
    TimedTask& LinearSystemFactorization()
    {
      return LinearSystemFactorization_;
    }
    TimedTask& LinearSystemBackSolve()
    {
      return LinearSystemBackSolve_;
    }
    TimedTask& LinearSystemStructureConverter()
    {
      return LinearSystemStructureConverter_;
    }
    TimedTask& LinearSystemStructureConverterInit()
    {
      return LinearSystemStructureConverterInit_;
    }
    TimedTask& QualityFunctionSearch()
    {
      return QualityFunctionSearch_;
    }
    TimedTask& TryCorrector()
    {
      return TryCorrector_;
    }

    TimedTask& Task1()
    {
      return Task1_;
    }
    TimedTask& Task2()
    {
      return Task2_;
    }
    TimedTask& Task3()
    {
      return Task3_;
    }
    TimedTask& Task4()
    {
      return Task4_;
    }
    TimedTask& Task5()
    {
      return Task5_;
    }
    TimedTask& Task6()
    {
      return Task6_;
    }
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
    /** Copy Constructor */
    TimingStatistics(const TimingStatistics&);

    /** Overloaded Equals Operator */
    void operator=(const TimingStatistics&);
    //@}

    /**@name All timed tasks. */
    //@{
    TimedTask OverallAlgorithm_;
    TimedTask PrintProblemStatistics_;
    TimedTask InitializeIterates_;
    TimedTask UpdateHessian_;
    TimedTask OutputIteration_;
    TimedTask UpdateBarrierParameter_;
    TimedTask ComputeSearchDirection_;
    TimedTask ComputeAcceptableTrialPoint_;
    TimedTask AcceptTrialPoint_;
    TimedTask CheckConvergence_;

    TimedTask PDSystemSolverTotal_;
    TimedTask PDSystemSolverSolveOnce_;
    TimedTask ComputeResiduals_;
    TimedTask LinearSystemScaling_;
    TimedTask LinearSystemSymbolicFactorization_;
    TimedTask LinearSystemFactorization_;
    TimedTask LinearSystemBackSolve_;
    TimedTask LinearSystemStructureConverter_;
    TimedTask LinearSystemStructureConverterInit_;
    TimedTask QualityFunctionSearch_;
    TimedTask TryCorrector_;

    TimedTask Task1_;
    TimedTask Task2_;
    TimedTask Task3_;
    TimedTask Task4_;
    TimedTask Task5_;
    TimedTask Task6_;
    //@}
  };

} // namespace Ipopt

#endif
