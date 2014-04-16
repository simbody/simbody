/*

Description: Provides a way to execute the same method on multiple instances
of a class in parallel on a thread pool with a syncronous
interface for invocation

Author: Nabeel Allana

*/

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <vector>
#include <deque>

namespace SimTK {

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
						_worker_thead_cond.wait_for(lock, std::chrono::duration<int, std::milli>(5));
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

			_tasks_remaining = _fnqueue.size();

			// wake up all workers
			_worker_thead_cond.notify_all();

			// don't return until all functions have been executed
			while (_tasks_remaining > 0)
			{
				_calling_thread_cond.wait_for(lock, 
					std::chrono::duration<int, std::milli>(5));
				continue;
			}

		}

		// Abort all threads and put class in invalid state
		void stop()
		{
			_sig_kill = true;

			for (std::size_t i = 0; i < _workers.size(); i++)
				_workers[i].join();

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

#endif