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

#include <iostream>
#include <string>
#include <algorithm>
#include <mutex>

using namespace std;

namespace SimTK {

static void threadBody(ThreadInfo& info);

ParallelExecutorImpl::ParallelExecutorImpl() : finished(false) {

    //By default, we use the total number of processors available of the
    //computer (including hyperthreads)
    numMaxThreads = ParallelExecutor::getNumProcessors();
    if(numMaxThreads <= 0)
      numMaxThreads = 1;
}
ParallelExecutorImpl::ParallelExecutorImpl(int numThreads) : finished(false) {

    // Set the maximum number of threads that we can use
    SimTK_APIARGCHECK_ALWAYS(numThreads > 0, "ParallelExecutorImpl",
                 "ParallelExecutorImpl", "Number of threads must be positive.");
    numMaxThreads = numThreads;
}
ParallelExecutorImpl::~ParallelExecutorImpl() {
    
    // Notify the threads that they should exit.
    
    std::unique_lock<std::mutex> lock(runMutex);
    finished = true;
    for (int i = 0; i < (int) threads.size(); ++i)
        threadInfo[i].running = true;
    runCondition.notify_all();
    lock.unlock();
    
    // Wait until all the threads have finished.
    
    for (int i = 0; i < (int) threads.size(); ++i)
        threads[i].join();
}
ParallelExecutorImpl* ParallelExecutorImpl::clone() const {
    return new ParallelExecutorImpl(numMaxThreads);
}
void ParallelExecutorImpl::execute(ParallelExecutor::Task& task, int times) {
  if (min(times, numMaxThreads) == 1) {
      //(1) NON-PARALLEL CASE:
      // Nothing is actually going to get done in parallel, so we might as well
      // just execute the task directly and save the threading overhead.
      task.initialize();
      for (int i = 0; i < times; ++i)
          task.execute(i);
      task.finish();
      return;
    }
    
    //(2) PARALLEL CASE:
    // We launch the maximum number of threads and save them for later use
    if(threads.size() < (size_t)numMaxThreads)
    {
      // We do not support numMaxThreads changing for a given instance of
      // ParallelExecutor.
      assert(threads.size() == 0);
      threadInfo.reserve(numMaxThreads);
      threads.resize(numMaxThreads);
      for (int i = 0; i < numMaxThreads; ++i) {
          threadInfo.emplace_back(i, this);
          threads[i] = std::thread(threadBody, std::ref(threadInfo[i]));
      }
    }
   
    // Initialize fields to execute the new task.
    std::unique_lock<std::mutex> lock(runMutex);
    currentTask = &task;
    currentTaskCount = times;
    waitingThreadCount = 0;
    for (int i = 0; i < (int)threadInfo.size(); ++i) {
        threadInfo[i].running = true;
    }

    // Wake up the worker threads and wait until they finish.
    runCondition.notify_all();
    waitCondition.wait(lock,
        [&] { return waitingThreadCount == (int)threads.size(); });
    lock.unlock();
}
void ParallelExecutorImpl::incrementWaitingThreads() {
    std::lock_guard<std::mutex> lock(runMutex);
    getCurrentTask().finish();
    waitingThreadCount++;
    if (waitingThreadCount == (int)threads.size()) {
        waitCondition.notify_one();
    }
}

thread_local bool ParallelExecutorImpl::isWorker(false);

/**
 * This function contains the code executed by the worker threads.
 */

void threadBody(ThreadInfo& info) {
    ParallelExecutorImpl::isWorker = true;
    ParallelExecutorImpl& executor = *info.executor;
    int threadCount = executor.getThreadCount();
    while (!executor.isFinished()) {
        
        // Wait for a Task to come in.
        
        std::unique_lock<std::mutex> lock(executor.getMutex());
        executor.getCondition().wait(lock, [&] { return info.running; });
        lock.unlock();
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
}

ParallelExecutor::ParallelExecutor() : HandleBase(new ParallelExecutorImpl()) {
}

ParallelExecutor::ParallelExecutor(int numThreads) : HandleBase(new ParallelExecutorImpl(numThreads)) {
}

ParallelExecutor* ParallelExecutor::clone() const{
    return new ParallelExecutor(getMaxThreads());
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
    int ncpu,retval;
    size_t len = 4;

    retval = sysctlbyname( "hw.physicalcpu", &ncpu, &len, NULL, 0 );
    if( retval == 0 ) {
       return (ncpu );
    } else {
       return(1);
    }
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
    return ParallelExecutorImpl::isWorker;
}
int ParallelExecutor::getMaxThreads() const{
    return getImpl().getMaxThreads();
}

} // namespace SimTK
