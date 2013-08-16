/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
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

#include "ParallelWorkQueueImpl.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include <pthread.h>
#include <utility>

using std::queue;

namespace SimTK {

static void* threadBody(void* args) {
    ParallelWorkQueueImpl& owner = *reinterpret_cast<ParallelWorkQueueImpl*>(args);
    queue<ParallelWorkQueue::Task*>& taskQueue = owner.updTaskQueue();
    pthread_mutex_t& queueLock = owner.getQueueLock();
    pthread_cond_t& waitForTaskCondition = owner.getWaitCondition();
    pthread_cond_t& queueFullCondition = owner.getQueueFullCondition();
    bool decrementTaskCount = false;
    while (!owner.isFinished() || !taskQueue.empty()) {
        pthread_mutex_lock(&queueLock);
        if (decrementTaskCount) {
            owner.markTaskCompleted();
            decrementTaskCount = false;
        }
        while (taskQueue.empty() && !owner.isFinished())
            pthread_cond_wait(&waitForTaskCondition, &queueLock);
        ParallelWorkQueue::Task* task = NULL;
        if (!taskQueue.empty()) {
            task = taskQueue.front();
            taskQueue.pop();
        }
        pthread_cond_signal(&queueFullCondition);
        pthread_mutex_unlock(&queueLock);
        if (task != NULL) {
            task->execute();
            delete task;
            decrementTaskCount = true;
        }
    }
    if (decrementTaskCount) {
        pthread_mutex_lock(&queueLock);
        owner.markTaskCompleted();
        pthread_mutex_unlock(&queueLock);
    }
    return 0;
}

ParallelWorkQueueImpl::ParallelWorkQueueImpl(int queueSize, int numThreads) : queueSize(queueSize), pendingTasks(0), finished(false) {
    pthread_mutex_init(&queueLock, NULL);
    pthread_cond_init(&waitForTaskCondition, NULL);
    pthread_cond_init(&queueFullCondition, NULL);
    threads.resize(numThreads);
    for (int i = 0; i < numThreads; ++i)
        pthread_create(&threads[i], NULL, threadBody, this);
}

ParallelWorkQueueImpl::~ParallelWorkQueueImpl() {
    // Wait for the theads to finish.

    pthread_mutex_lock(&queueLock);
    finished = true;
    pthread_cond_broadcast(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
    for (int i = 0; i < (int) threads.size(); ++i)
        pthread_join(threads[i], NULL);

    // Clean up memory.
    
    pthread_mutex_destroy(&queueLock);
    pthread_cond_destroy(&waitForTaskCondition);
    pthread_cond_destroy(&queueFullCondition);
}

ParallelWorkQueueImpl* ParallelWorkQueueImpl::clone() const {
    return new ParallelWorkQueueImpl(queueSize, threads.size());
}

void ParallelWorkQueueImpl::addTask(ParallelWorkQueue::Task* task) {
    pthread_mutex_lock(&queueLock);
    while ((int)taskQueue.size() >= queueSize)
        pthread_cond_wait(&queueFullCondition, &queueLock);
    taskQueue.push(task);
    ++pendingTasks;
    pthread_cond_signal(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
}

void ParallelWorkQueueImpl::flush() {
    pthread_mutex_lock(&queueLock);
    while (pendingTasks > 0)
       pthread_cond_wait(&queueFullCondition, &queueLock);
    pthread_mutex_unlock(&queueLock);
}

queue<ParallelWorkQueue::Task*>& ParallelWorkQueueImpl::updTaskQueue() {
    return taskQueue;
}

bool ParallelWorkQueueImpl::isFinished() const {
    return finished;
}

pthread_mutex_t& ParallelWorkQueueImpl::getQueueLock() {
    return queueLock;
}

pthread_cond_t& ParallelWorkQueueImpl::getWaitCondition() {
    return waitForTaskCondition;
}

pthread_cond_t& ParallelWorkQueueImpl::getQueueFullCondition() {
    return queueFullCondition;
}

void ParallelWorkQueueImpl::markTaskCompleted() {
    --pendingTasks;
    pthread_cond_signal(&queueFullCondition);
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
