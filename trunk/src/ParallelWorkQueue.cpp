/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
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
    while (!owner.isFinished() || !taskQueue.empty()) {
        pthread_mutex_lock(&queueLock);
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
        }
    }
    return 0;
}

ParallelWorkQueueImpl::ParallelWorkQueueImpl(int queueSize, int numThreads) : queueSize(queueSize), finished(false) {
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
    SimTK_ASSERT_ALWAYS(!finished, "Tried to add a Task after calling finish()");
    pthread_mutex_lock(&queueLock);
    while ((int)taskQueue.size() >= queueSize)
        pthread_cond_wait(&queueFullCondition, &queueLock);
    taskQueue.push(task);
    pthread_cond_signal(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
}

void ParallelWorkQueueImpl::flush() {
    pthread_mutex_lock(&queueLock);
    while (taskQueue.size() > 0)
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

ParallelWorkQueue::ParallelWorkQueue(int queueSize, int numThreads) : HandleBase(new ParallelWorkQueueImpl(queueSize, numThreads)) {
}

void ParallelWorkQueue::addTask(Task* task) {
    updImpl().addTask(task);
}

void ParallelWorkQueue::flush() {
    updImpl().flush();
}

} // namespace SimTK
