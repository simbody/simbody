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

#include "ParallelExecutorImpl.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include <pthread.h>

namespace SimTK {

static void* threadBody(void* args);

ParallelExecutorImpl::ParallelExecutorImpl(int numThreads) : finished(false) {

    // Construct all the threading related objects we will need.
    
    SimTK_APIARGCHECK_ALWAYS(numThreads > 0, "ParallelExecutorImpl", "ParallelExecutorImpl", "Number of threads must be positive.");

    threads.resize(numThreads);
    pthread_mutex_init(&runLock, NULL);
    pthread_cond_init(&runCondition, NULL);
    pthread_cond_init(&waitCondition, NULL);
    for (int i = 0; i < numThreads; ++i) {
        threadInfo.push_back(new ThreadInfo(i, this));
        pthread_create(&threads[i], NULL, threadBody, threadInfo[i]);
    }
}
ParallelExecutorImpl::~ParallelExecutorImpl() {
    
    // Notify the threads that they should exit.
    
    pthread_mutex_lock(&runLock);
    finished = true;
    for (int i = 0; i < (int) threads.size(); ++i)
        threadInfo[i]->running = true;
    pthread_cond_broadcast(&runCondition);
    pthread_mutex_unlock(&runLock);
    
    // Wait until all the threads have finished.
    
    for (int i = 0; i < (int) threads.size(); ++i)
        pthread_join(threads[i], NULL);
    
    // Clean up threading related objects.
    
    pthread_mutex_destroy(&runLock);
    pthread_cond_destroy(&runCondition);
    pthread_cond_destroy(&waitCondition);
}
ParallelExecutorImpl* ParallelExecutorImpl::clone() const {
    return new ParallelExecutorImpl(threads.size());
}
void ParallelExecutorImpl::execute(ParallelExecutor::Task& task, int times) {
    if (times == 1 || threads.size() == 1) {
        // Nothing is actually going to get done in parallel, so we might as well
        // just execute the task directly and save the threading overhead.
        
        task.initialize();
        for (int i = 0; i < times; ++i)
            task.execute(i);
        task.finish();
        return;
    }
    
    // Initialize fields to execute the new task.
    
    pthread_mutex_lock(&runLock);
    currentTask = &task;
    currentTaskCount = times;
    waitingThreadCount = 0;
    for (int i = 0; i < (int) threadInfo.size(); ++i)
        threadInfo[i]->running = true;

    // Wake up the worker threads and wait until they finish.
    
    pthread_cond_broadcast(&runCondition);
    do {
        pthread_cond_wait(&waitCondition, &runLock);
    } while (waitingThreadCount < (int) threads.size());
    pthread_mutex_unlock(&runLock);
}
void ParallelExecutorImpl::incrementWaitingThreads() {
    pthread_mutex_lock(&runLock);
    getCurrentTask().finish();
    waitingThreadCount++;
    if (waitingThreadCount == (int)threads.size()) {
        pthread_cond_signal(&waitCondition);
    }
    pthread_mutex_unlock(&runLock);
}

ThreadLocal<bool> ParallelExecutorImpl::isWorker(false);

/**
 * This function contains the code executed by the worker threads.
 */

void* threadBody(void* args) {
    ParallelExecutorImpl::isWorker.upd() = true;
    ThreadInfo& info = *reinterpret_cast<ThreadInfo*>(args);
    ParallelExecutorImpl& executor = *info.executor;
    int threadCount = executor.getThreadCount();
    while (!executor.isFinished()) {
        
        // Wait for a Task to come in.
        
        pthread_mutex_lock(executor.getLock());
        while (!info.running) {
            pthread_cond_wait(executor.getCondition(), executor.getLock());
        }
        pthread_mutex_unlock(executor.getLock());
        if (!executor.isFinished()) {
            
            // Execute the task for all the indices belonging to this thread.
            
            int count = executor.getCurrentTaskCount();
            ParallelExecutor::Task& task = executor.getCurrentTask();
            task.initialize();
            int index = info.index;
            try {
                while (index < count) {
                    task.execute(index);
                    index += threadCount;
                }
            }
            catch (const std::exception& ex) {
                std::cerr <<"The parallel task threw an unhandled exception:"<< std::endl;
                std::cerr <<ex.what()<< std::endl;
            }
            catch (...) {
                std::cerr <<"The parallel task threw an error."<< std::endl;
            }
            info.running = false;
            executor.incrementWaitingThreads();
        }
    }
    delete &info;
    return 0;
}

ParallelExecutor::ParallelExecutor(int numThreads) : HandleBase(new ParallelExecutorImpl(numThreads)) {
}

void ParallelExecutor::execute(Task& task, int times) {
    updImpl().execute(task, times);
}

#ifdef __APPLE__
   #include <sys/sysctl.h>
   #include <dlfcn.h>
#elif _WIN32
   #include <windows.h>
#elif __linux__
   #include <dlfcn.h>
   #include <unistd.h>
#else
  #error "Architecture unsupported"
#endif

int ParallelExecutor::getNumProcessors() {
#ifdef __APPLE__
    int numCPU = 1; //by default set the number of CPUs to 1
    int mib[4];
    size_t len = sizeof(numCPU); 

    /* set the mib for hw.ncpu */
    mib[0] = CTL_HW;
    mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;

    /* get the number of CPUs from the system */
    sysctl(mib, 2, &numCPU, &len, NULL, 0);

    if( numCPU < 1 ) 
    {
         mib[1] = HW_NCPU;
         sysctl( mib, 2, &numCPU, &len, NULL, 0 );

         if( numCPU < 1 )
         {
              numCPU = 1;
         }
    }
    return numCPU;
#else
#ifdef __linux__
    long nProcessorsOnline     = sysconf(_SC_NPROCESSORS_ONLN);
    if( nProcessorsOnline == -1 )  {
        return(1);
    } else {
        return( (int)nProcessorsOnline );
    }
#else
    // Windows
    SYSTEM_INFO siSysInfo;
    int ncpu;

    // Copy the hardware information to the SYSTEM_INFO structure.

    GetSystemInfo(&siSysInfo);

    // Display the contents of the SYSTEM_INFO structure.

    ncpu =  siSysInfo.dwNumberOfProcessors;
    if( ncpu < 1 ) ncpu = 1;
    return(ncpu);
#endif
#endif
}

bool ParallelExecutor::isWorkerThread() {
    return ParallelExecutorImpl::isWorker.get();
}

} // namespace SimTK
