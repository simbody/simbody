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

#include "ParallelWorkQueueImpl.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include <utility>
#include <mutex>
#include <condition_variable>

using std::queue;

namespace SimTK {

static void threadBody(ParallelWorkQueueImpl& owner) {
    queue<ParallelWorkQueue::Task*>& taskQueue = owner.updTaskQueue();
    std::mutex& queueMutex = owner.getQueueMutex();
    std::condition_variable& waitForTaskCondition = owner.getWaitCondition();
    std::condition_variable& queueFullCondition = owner.getQueueFullCondition();
    bool decrementTaskCount = false;
    while (!owner.isFinished() || !taskQueue.empty()) {
        std::unique_lock<std::mutex> lock(queueMutex);
        if (decrementTaskCount) {
            owner.markTaskCompleted();
            decrementTaskCount = false;
        }
        waitForTaskCondition.wait(lock,
                [&] { return !taskQueue.empty() || owner.isFinished(); });
        ParallelWorkQueue::Task* task = NULL;
        if (!taskQueue.empty()) {
            task = taskQueue.front();
            taskQueue.pop();
        }
        queueFullCondition.notify_one();
        lock.unlock();
        if (task != NULL) {
            task->execute();
            delete task;
            decrementTaskCount = true;
        }
    }
    if (decrementTaskCount) {
        std::lock_guard<std::mutex> lock(queueMutex);
        owner.markTaskCompleted();
    }
}

ParallelWorkQueueImpl::ParallelWorkQueueImpl(int queueSize, int numThreads) : queueSize(queueSize), pendingTasks(0), finished(false) {
    threads.resize(numThreads);
    for (int i = 0; i < numThreads; ++i)
        threads[i] = std::thread(threadBody, std::ref(*this));
}

ParallelWorkQueueImpl::~ParallelWorkQueueImpl() {
    // Wait for the theads to finish.

    std::unique_lock<std::mutex> lock(queueMutex);
    finished = true;
    waitForTaskCondition.notify_all();
    lock.unlock();
    for (int i = 0; i < (int) threads.size(); ++i)
        threads[i].join();
}

ParallelWorkQueueImpl* ParallelWorkQueueImpl::clone() const {
    return new ParallelWorkQueueImpl(queueSize, threads.size());
}

void ParallelWorkQueueImpl::addTask(ParallelWorkQueue::Task* task) {
    std::unique_lock<std::mutex> lock(queueMutex);
    queueFullCondition.wait(lock,
            [this] { return (int)taskQueue.size() < queueSize; });
    taskQueue.push(task);
    ++pendingTasks;
    waitForTaskCondition.notify_one();
    lock.unlock();
}

void ParallelWorkQueueImpl::flush() {
    std::unique_lock<std::mutex> lock(queueMutex);
    queueFullCondition.wait(lock, [this] { return pendingTasks == 0; });
    lock.unlock();
}

queue<ParallelWorkQueue::Task*>& ParallelWorkQueueImpl::updTaskQueue() {
    return taskQueue;
}

bool ParallelWorkQueueImpl::isFinished() const {
    return finished;
}

std::mutex& ParallelWorkQueueImpl::getQueueMutex() {
    return queueMutex;
}

std::condition_variable& ParallelWorkQueueImpl::getWaitCondition() {
    return waitForTaskCondition;
}

std::condition_variable& ParallelWorkQueueImpl::getQueueFullCondition() {
    return queueFullCondition;
}

void ParallelWorkQueueImpl::markTaskCompleted() {
    --pendingTasks;
    queueFullCondition.notify_one();
}

ParallelWorkQueue::ParallelWorkQueue(int queueSize, int numThreads) : HandleBase(new ParallelWorkQueueImpl(queueSize, numThreads)) {
}

void ParallelWorkQueue::addTask(Task* task) {
    updImpl().addTask(task);
}

void ParallelWorkQueue::flush() {
    updImpl().flush();
}

} // namespace SimTK
