#ifndef SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE11_H_
#define SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE11_H_

/* -------------------------------------------------------------------------- *
*                       Simbody(tm): SimTKcommon                             *
* -------------------------------------------------------------------------- *
* This is part of the SimTK biosimulation toolkit originating from           *
* Simbios, the NIH National Center for Physics-Based Simulation of           *
* Biological Structures at Stanford, funded under the NIH Roadmap for        *
* Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
*                                                                            *
* Portions copyright (c) 2014 Stanford University and the Authors.           *
* Authors: Nabeel Allana                                                     *
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

#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <vector>
#include <deque>

namespace SimTK {

/**
* This class is used for performing multithreaded computations.\ It maintains a queue of tasks to be
* executed, and a pool of threads for executing them.  This class is similar to the ParallelWorkQueue class
* with the difference being that this class uses C++11 threading, whereas the ParallelWorkQueue uses
* the traditional p-threads libaray. This class is, as a result, more flexible than the ParallelWorkQueue
* class in its usage. Usage:
*
* <pre>
* void foo() { doStuff(); }
* void foo2(int a) { compute(a); }
*
* class bar { void exec() { doMoreStuff(); } };
*
* ParallelWorkQueue11 queue(20);
* bar myBar; 
*
* // plain old functions
* queue.push(&foo);
*
* // functions with pre-bound arguments
* queue.push(std::bind(&foo, 100));
*
* // member functions
* queue.push(std::bind(&bar::exec, &myBar));
*
* // execute everything in queue. returns when all are complete
* queue.execute();
*
* </pre>
*
* Each function will be called by the queue. A few examples on loading methods in the queue are provided.
* For more, see std::function and std::bind in the standard library. The invocations are
* done in parallel on multiple threads, so you cannot make any assumptions about what order they will occur in
* or which ones will happen at the same time.
*
* The threads are created in the ParallelWorkQueue's constructor and remain active until it is deleted.
* This means that creating a ParallelWorkQueue is a somewhat expensive operation, but it may then be
* used repeatedly for executing various calculations. 
*/
class ParallelWorkQueue11
{
private:
	std::deque<std::function<void()>> _fnqueue;

	std::mutex _execute_mutex;
	std::mutex _queue_mutex;
	std::mutex _task_len_mutex;
	std::mutex _calling_thread_mtx;

	std::vector<std::thread> _workers;

	std::condition_variable _worker_thead_cond;
	std::condition_variable _calling_thread_cond;

	bool _sig_kill;
	bool _sync_return;
	int _tasks_remaining;

	// Internal worker method - each thread runs an
	// instance of this method
	void worker()
	{
		std::function<void()> fn;
		while (!_sig_kill)
		{
			{
				std::unique_lock<std::mutex> lock(_queue_mutex);

				if (_fnqueue.empty())
				{
					_worker_thead_cond.wait(lock);
					continue;
				}

				fn = _fnqueue.front();
				_fnqueue.pop_front();
			}
			try
			{
				fn();
			}
			catch (std::exception& ex)
			{
				_sig_kill = true;
				throw ex;
			}
			if (!_sig_kill)
			{
				std::unique_lock<std::mutex> lock(_task_len_mutex);

				--_tasks_remaining;

				_calling_thread_cond.notify_one();
			}
		}

	}
public:

	void push(std::function<void()> fn)
	{
		if (_sig_kill)
			throw std::exception("Thread pool has been stopped.");

		std::unique_lock<std::mutex> exec(_execute_mutex);
		std::unique_lock<std::mutex> lock(_queue_mutex);

		_fnqueue.push_back(fn);
	}

	void execute()
	{
		if (_sig_kill)
			throw std::exception("Thread pool has been stopped.");

		// prevent anyone else from calling execute
		// until queue is flushed.
		std::unique_lock<std::mutex>
			exec(_execute_mutex);
		std::unique_lock<std::mutex>
			lock(_calling_thread_mtx);
		std::unique_lock<std::mutex>
			tlLock(_task_len_mutex);

		_tasks_remaining = _fnqueue.size();

		// wake up all workers
		_worker_thead_cond.notify_all();

		// don't return until all functions have been executed
		while (_tasks_remaining > 0)
		{
			_calling_thread_cond.wait_for(tlLock,
				std::chrono::duration<int, std::milli>(5));
			continue;
		}

	}

	// Change the number of worker threads
	void resize(std::size_t n)
	{
		stop();
		_sig_kill = false;

		for (std::size_t i = 0; i < n; i++)
			_workers.push_back(std::thread(std::bind(&ParallelWorkQueue11::worker, this)));
	}

	// Abort all threads and put class in invalid state
	void stop()
	{
		_sig_kill = true;

		_worker_thead_cond.notify_all();

		for (std::size_t i = 0; i < _workers.size(); i++)
			_workers[i].join();

		_workers.clear();

	}

	// Initialize and spawn worker threads
	ParallelWorkQueue11(std::size_t nThreads, bool async = false)
		: _sig_kill(false),
		_sync_return(!async)
	{
		for (std::size_t i = 0; i < nThreads; i++)
			_workers.push_back(std::thread(std::bind(&ParallelWorkQueue11::worker, this)));
	}

	~ParallelWorkQueue11()
	{
		stop();
	}
};

}

#endif // SimTK_SimTKCOMMON_PARALLEL_WORK_QUEUE11_H_