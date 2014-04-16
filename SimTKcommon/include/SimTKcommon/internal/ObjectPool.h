
#include <deque>
#include <mutex>
#include <condition_variable>
#include <list>

template<class TObj>
class ObjectPool
{
public:
	~ObjectPool()
	{
		while (!_objs.empty())
		{
			delete _objs.front();
			_objs.pop_front();
		}
	}
	template<class...params>
	void initialize(int nInstances, params... p)
	{
		for (int i = 0; i < nInstances; ++i)
		{
			TObj* ob = new TObj(p...);

			_obl.push_back(ob);
			_objs.push_back(ob);
		}
	}
	TObj* get()
	{
		TObj* rtn;
		std::unique_lock<std::mutex>
			lock(_objs_mtx);
		
		while (_objs.empty())
			_cond.wait(lock);

		rtn = _objs.front();
		_objs.pop_front();

		return rtn;

	}
	const std::list<TObj*>& list() const
	{
		return _obl;
	}
	void push(TObj* o)
	{
		std::unique_lock<std::mutex>
			lock(_objs_mtx);
		_objs.push_back(o);
		_cond.notify_one();
	}
private:
	std::list<TObj*> _obl;
	std::deque<TObj*> _objs;
	std::condition_variable _cond;
	std::mutex _objs_mtx;
};