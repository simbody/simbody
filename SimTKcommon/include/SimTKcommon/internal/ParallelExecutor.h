#ifndef SimTK_SimTKCOMMON_PARALLEL_EXECUTOR_H_
#define SimTK_SimTKCOMMON_PARALLEL_EXECUTOR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "PrivateImplementation.h"
#include <thread>
#include <iostream>
 
namespace SimTK {

class ParallelExecutor;
class ParallelExecutorImpl;

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the ParallelExecutor class but should not be defined any other time.

#ifndef SimTK_SIMTKCOMMON_DEFINING_PARALLEL_EXECUTOR
    extern template class PIMPLHandle<ParallelExecutor, ParallelExecutorImpl>;
#endif

/**
 * This class is used for performing multithreaded computations.  To use it,
 * define a subclass of ParallelExecutor::Task that performs some computation.
 * Then create a ParallelExecutor object and ask it to execute the task:
 * 
 * <pre>
 * ParallelExecutor executor;
 * executor.execute(myTask, times);
 * </pre>
 * 
 * The Task's execute() method will be called the specified number of times,
 * with each invocation being given a different index value from 0 to times-1.
 * The invocations are done in parallel on multiple threads, so you cannot make
 * any assumptions about what order they will occur in or which ones will
 * happen at the same time.
 * 
 * The threads are created in the ParallelExecutor's constructor and remain
 * active until it is deleted. This means that creating a ParallelExecutor is a
 * somewhat expensive operation, but it may then be used repeatedly for
 * executing various calculations.  By default, the number of threads is chosen
 * to be equal to the number of available processor cores.  You can optionally
 * specify a different number of threads to create.  For example, using more
 * threads than processors can sometimes lead to better processor
 * utilitization.  Alternatively, if the Task will only be executed four times,
 * you might specify min(4, ParallelExecutor::getNumProcessors()) to avoid
 * creating extra threads that will never have any work to do.
 *
 * You may find it useful to use "thread local" variables with your parallel
 * tasks. A thread local variable may have a different value on each thread
 * Thread-local storage is useful in many situations when writing multithreaded
 * code. For example, it can be used as temporary workspace for calculations.
 * If a single workspace object were created, all access to it would need to be
 * synchronized to prevent threads from overwriting each other's values. Making
 * a variable `thread_local T` means that a separate workspace object of type
 * `T` will automatically be created for each thread. That object will have
 * "thread scope" meaning it will be destructed only on thread termination.
 * 
 * One way to use thread local variables is to declare a static thread local
 * member variable in your ParallelExecutor::Task:
 * @code{.cpp}
 * class MyTask : public ParallelExecutor::Task {
 *     ...
 * public:
 *     void execute(int index) override {
 *         workspace += 1.0; // incremented separately on each thread.
 *     }
 * private:
 *     static thread_local double workspace;
 * };
 * // Initialize the static variable out-of-line.
 * thread_local double workspace = 0;
 * @endcode
 */

class SimTK_SimTKCOMMON_EXPORT ParallelExecutor : public PIMPLHandle<ParallelExecutor, ParallelExecutorImpl> {
public:
    class Task;
    /**
     * Construct a ParallelExecutor. By default, constructs a ParallelExecutor
     * with the number of threads equal to the number of total processors on the
     * computer (including hyperthreads).
     */
    explicit ParallelExecutor();
    /**
     * Construct a ParallelExecutor.
     * 
     * @param maxThreads the maximum number of threads that the ParallelExecutor
     * is allowed to launch
     */
    explicit ParallelExecutor(int maxThreads);
    /**
     * Clone the ParallelExecutor.
     *
     * @return Returns a ParallelExecutor with the same number of maxThreads
     */
    ParallelExecutor* clone() const;
    /**
     * Execute a parallel task.
     * 
     * @param task    the Task to execute
     * @param times   the number of times the Task should be executed
     */
    void execute(Task& task, int times);
    /**
     * Get the total number of available processor cores (physical cores and
     * hyperthreads on Intel architecture). If the number of threads is not
     * computable or well defined, the function returns 0.
     */
    static int getNumProcessors();
    /**
     * Determine whether the thread invoking this method is a worker thread
     * created by ParallelExecutor.
     */
    static bool isWorkerThread();
    /**
     * Get the maximum number of thread contexts that the ParallelExecutor is
     * currently allowed to use.
     */
    int getMaxThreads() const;
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
