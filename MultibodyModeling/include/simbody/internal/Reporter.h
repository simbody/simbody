#ifndef SIMTK_SIMBODY_REPORTER_H_
#define SIMTK_SIMBODY_REPORTER_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * User-visible, client side declaration of Reporter.
 * Concrete objects derived from Reporter are able to capture results from
 * "solvers", where that term should be understood as any producer of States.
 * The (Model,State) pair is passed to the Reporter's report() method; after
 * that it's up to you ...
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/State.h"
#include "simbody/internal/Model.h"

namespace simtk {

class Reporter {
public:
    Reporter();
    ~Reporter();
    Reporter(const Reporter&);
    Reporter& operator=(const Reporter&);

    enum Reason {
        BeforeStudy,
        AfterStudy,
        StepIntervalReached,
        TimeIntervalReached,
        AdHocReport
    };

    void report(const Model&, const State&, Reason=AdHocReport) const;

    static String getReasonName(Reason);

private:
    class ReporterRep* rep;
    friend class ReporterRep;
};

/**
 * This class specifies when we would like our Reporter to be called.
 * XXX Later this should be unified with the Event handling system.
 */
class ReportInterval {
    enum {
        BeforeStudy = 0x001,
        AfterStudy  = 0x002,
        StepInterval = 0x004,
        TimeInterval = 0x008
    };

public:
    /// Default constructor yields no reporting at all; use add() methods below.
    ReportInterval();

    /// This is the commonly used constructor which provides before/after reporting
    /// and time interval and/or step interval reporting (set dt or ds to 0 to
    /// suppress time or step reporting, resp.).
    explicit ReportInterval(const Real& dt, unsigned int ds=0); 

         
    void addBeforeStudy() { when |= BeforeStudy; }
    void addAfterStudy()  { when |= AfterStudy; }
    void addStepInterval(int dstep, bool measureFromZero=false)
        { when |= StepInterval; stepInterval=dstep; stepFromZero=measureFromZero; }
    void addTimeInterval(int dt, bool measureFromZero=false)
        { when |= TimeInterval; timeInterval=dt; timeFromZero=measureFromZero; }

    bool isReportingStep(unsigned int startStep, unsigned int curStep) const;
    bool isReportingTime(const Real& startTime,
                         const Real& lastReportTime, 
                         const Real& curTime) const;
    Real calcNextReportingTime(const Real& startTime, const Real& lastReportTime) const;
    
    bool shouldReportBeforeStudy() const { return (when & BeforeStudy)!=0; }
    bool shouldReportAfterStudy()  const { return (when & AfterStudy)!=0; }
    bool shouldReportAtTimeInterval() const { return (when & TimeInterval)!=0; }
    bool shouldReportAtStepInterval() const { return (when & StepInterval)!=0; }
    bool isTimeMeasuredFromZero() const { return timeFromZero; }
    bool isStepMeasuredFromZero() const { return stepFromZero; }
    unsigned int getStepInterval() const { return stepInterval; }
    Real         getTimeInterval() const { return timeInterval; }  
  
private:
    unsigned int when;
    bool         stepFromZero;  // i.e., report when (step%interval)==0 (or time)
    bool         timeFromZero;  //   default is ((step-startStep)%interval)==0
    unsigned int stepInterval;
    Real         timeInterval;
};

} // namespace simtk

#endif // SIMTK_SIMBODY_REPORTER_H_
