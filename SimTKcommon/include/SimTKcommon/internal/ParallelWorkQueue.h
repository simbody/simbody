#ifndef SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_H_
#define SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "ParallelExecutor.h"
#include "PrivateImplementation.h"

namespace SimTK {

class ParallelWorkQueue;
class ParallelWorkQueueImpl;

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the ParallelWorkQueue class but should not be defined any other time.
#ifndef SimTK_SIMTKCOMMON_DEFINING_PARALLEL_WORK_QUEUE
    extern template class PIMPLHandle<ParallelWorkQueue, ParallelWorkQueueImpl>;
#endif

/**
 * This class is used for performing multithreaded computations.\ It maintains a queue of tasks to be
 * executed, and a pool of threads for executing them.  To use it, define one or more subclasses of
 * ParallelWorkQueue::Task that performs computations.  Then create a ParallelWorkQueue and add Tasks
 * to it:
 *
 * <pre>
 * ParallelWorkQueue queue(20);
 * queue.addTask(new MyTask());
 * </pre>
 *
 * Each Task's execute() method will be called, and then the Task will be deleted.  The invocations are
 * done in parallel on multiple threads, so you cannot make any assumptions about what order they will occur in
 * or which ones will happen at the same time.
 *
 * The threads are created in the ParallelWorkQueue's constructor and remain active until it is deleted.
 * This means that creating a ParallelWorkQueue is a somewhat expensive operation, but it may then be
 * used repeatedly for executing various calculations.  By default, the number of threads is chosen
 * to be equal to the number of available processor cores.  You can optionally specify a different number
 * of threads to create.  For example, using more threads than processors can sometimes lead to better
 * processor utilitization.  Alternatively, if only four Tasks will be executed, you might specify
 * min(4, ParallelExecutor::getNumProcessors()) to avoid creating extra threads that will never
 * have any work to do.
 */

class SimTK_SimTKCOMMON_EXPORT ParallelWorkQueue : public PIMPLHandle<ParallelWorkQueue, ParallelWorkQueueImpl> {
public:
    class Task;
    /**
     * Construct a ParallelWorkQueue.
     *
     * @param queueSize  the maximum number of Tasks that can be in the queue waiting to start executing at any time
     * @param numThreads the number of threads to create.  By default, this is set equal to the number of processors.
     */
    explicit ParallelWorkQueue(int queueSize, int numThreads = ParallelExecutor::getNumProcessors());
    /**
     * Add a Task to the queue.  If the queue is currently full, this method will block until space is freed up in
     * the queue by the worker threads.  The queue assumes ownership of the Task object and deletes it once it has
     * finished executing.
     *
     * @param task      the Task to add to the queue
     */
    void addTask(Task* task);
    /**
     * Block until all Tasks that have been added to the queue finish executing.
     */
    void flush();
};

/**
 * Concrete subclasses of this abstract class represent tasks that can be executed by a ParallelWorkQueue.
 */

class ParallelWorkQueue::Task {
public:
    virtual ~Task() {
    }
    /**
     * This method defines the task to be performed.  At some point after the Task is added to a
     * ParallelWorkQueue, this method will be called by one of the worker threads.
     */
    virtual void execute() = 0;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_H_
