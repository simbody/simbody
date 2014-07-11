#ifndef SimTK_SimTKCOMMON_OBJECT_POOL_H_
#define SimTK_SimTKCOMMON_OBJECT_POOL_H_

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

#include <deque>
#include <mutex>
#include <condition_variable>
#include <list>

namespace SimTK
{

template<class TObj>
class ObjectPool
{
public:
	~ObjectPool()
	{
		destroy();
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
	template<class...params>
	void resize(int nInstances, params... p)
	{
		destroy();
		initialize<params>(nInstances, p...);
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

	void destroy()
	{
		while (!_objs.empty())
		{
			delete _objs.front();
			_objs.pop_front();
		}
	}

	std::list<TObj*> _obl;
	std::deque<TObj*> _objs;
	std::condition_variable _cond;
	std::mutex _objs_mtx;
};

}

#endif // SimTK_SimTKCOMMON_OBJECT_POOL_H_