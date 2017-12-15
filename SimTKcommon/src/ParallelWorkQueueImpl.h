#ifndef SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_IMPL_H_
#define SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-17 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Christopher Dembia                                           *
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

#include "SimTKcommon/internal/ParallelWorkQueue.h"
#include "SimTKcommon/internal/Array.h"
#include <queue>
#include <mutex>
#include <condition_variable>

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
    std::mutex& getQueueMutex();
    std::condition_variable& getWaitCondition();
    std::condition_variable& getQueueFullCondition();
    void markTaskCompleted();
private:
    const int queueSize;
    int pendingTasks;
    bool finished;
    std::queue<ParallelWorkQueue::Task*> taskQueue;
    std::mutex queueMutex;
    std::condition_variable waitForTaskCondition, queueFullCondition;
    SimTK::Array_<std::thread> threads;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE_IMPL_H_
