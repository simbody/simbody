#ifndef SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_IMPL_H_
#define SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_IMPL_H_

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

#include "SimTKcommon/internal/ParallelWorkQueue.h"
#include "SimTKcommon/internal/Array.h"
#include <pthread.h>
#include <queue>

namespace SimTK {

class ParallelExecutor;

/**
 * This is the internal implementation class for ParallelWorkQueue.
 */

class ParallelWorkQueueImpl : public PIMPLImplementation<ParallelWorkQueue, ParallelWorkQueueImpl> {
public:
    class ExecutorTask;
    ParallelWorkQueueImpl(int queueSize, int numThreads);
    ~ParallelWorkQueueImpl();
    ParallelWorkQueueImpl* clone() const;
    void addTask(ParallelWorkQueue::Task* task);
    void flush();
    std::queue<ParallelWorkQueue::Task*>& updTaskQueue();
    bool isFinished() const;
    pthread_mutex_t& getQueueLock();
    pthread_cond_t& getWaitCondition();
    pthread_cond_t& getQueueFullCondition();
    void markTaskCompleted();
private:
    const int queueSize;
    int pendingTasks;
    bool finished;
    std::queue<ParallelWorkQueue::Task*> taskQueue;
    pthread_mutex_t queueLock;
    pthread_cond_t waitForTaskCondition, queueFullCondition;
    Array_<pthread_t> threads;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_IMPL_H_
