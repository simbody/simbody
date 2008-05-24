#ifndef SimTK_SimTKCOMMON_PARALLEL_EXECUTOR_H_
#define SimTK_SimTKCOMMON_PARALLEL_EXECUTOR_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "PrivateImplementation.h"

namespace SimTK {

class ParallelExecutor;
class ParallelExecutorImpl;

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the ParallelExecutor class but should not be defined any other time.
#ifndef SimTK_SIMTKCOMMON_DEFINING_PARALLEL_EXECUTOR
    extern template class PIMPLHandle<ParallelExecutor, ParallelExecutorImpl>;
#endif

/**
 * This class is used for performing multithreaded computations.  To use it, define a subclass of
 * ParallelExecutor::Task that performs some computation.  Then create a ParallelExecutor object
 * and ask it to execute the task:
 * 
 * <pre>
 * ParallelExecutor executor;
 * executor.execute(myTask, times);
 * </pre>
 * 
 * The Task's execute() method will be called the specified number of times, with each invocation
 * being given a different index value from 0 to times-1.  The invocations are done in parallel
 * on multiple threads, so you cannot make any assumptions about what order they will occur in
 * or which ones will happen at the same time.
 * 
 * The threads are created in the ParallelExecutor's constructor and remain active until it is deleted.
 * This means that creating a ParallelExecutor is a somewhat expensive operation, but it may then be
 * used repeatedly for executing various calculations.  By default, the number of threads is chosen
 * to be equal to the number of available processor cores.  You can optionally specify a different number
 * of threads to create.  For example, using more threads than processors can sometimes lead to better
 * processor utilitization.  Alternatively, if the Task will only be executed four times, you might
 * specify max(4, ParallelExecutor::getNumProcessors()) to avoid creating extra threads that will never
 * have any work to do.
 */

class SimTK_SimTKCOMMON_EXPORT ParallelExecutor : public PIMPLHandle<ParallelExecutor, ParallelExecutorImpl> {
public:
    class Task;
    /**
     * Construct a ParallelExecutor.
     * 
     * @param numThreads the number of threads to create.  By default, this is set equal to the number
     * of processors.
     */
    ParallelExecutor(int numThreads = getNumProcessors());
    /**
     * Execute a parallel task.
     * 
     * @param task    the Task to execute
     * @param times   the number of times the Task should be executed
     */
    void execute(Task& task, int times);
    /**
     * Get the number of available processor cores.
     */
    static int getNumProcessors();
};

/**
 * Concrete subclasses of this abstract class represent tasks that can be executed by a ParallelExecutor.
 */

class ParallelExecutor::Task {
public:
    virtual ~Task() {
    }
    /**
     * This method defines the task to be performed.  When the Task is passed to a ParallelExecutor's execute()
     * method, this method will be called in parallel the specified number of times.
     */
    virtual void execute(int index) = 0;
    /**
     * This method is invoked once by each worker thread before the task is executed.  This can be used to
     * initialize thread-local storage.
     */
    virtual void initialize() {
    }
    /**
     * This method is invoked once by each worker thread after all invocations of the task on that thread are complete.
     * This can be used to clean up thread-local storage, or to record per-thread results.  All calls to this method
     * are synchronized, so it can safely write to global variables without danger of interference between worker threads.
     */
    virtual void finish() {
    }
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_EXECUTOR_H_
