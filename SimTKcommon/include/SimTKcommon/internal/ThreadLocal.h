#ifndef SimTK_SimTKCOMMON_THREAD_LOCAL_H_
#define SimTK_SimTKCOMMON_THREAD_LOCAL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-15 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include <pthread.h>
#include <map>
#include <set>

namespace SimTK {

/** <b>(Deprecated)</b> This class represents a "thread local" variable: one 
which may have a different value on each thread; use C++11 `thread_local`
instead. This class is no longer necessary since C++11 has `thread_local` as a 
built-in keyword -- you should use that instead.

Thread-local storage is useful in many situations when writing multithreaded 
code. For example, it can be used as temporary workspace for calculations. If a 
single workspace object were created, all access to it would need to be 
synchronized to prevent threads from overwriting each other's values. Using a 
`ThreadLocal<T>` instead means that a separate workspace object of type
`T` will automatically be created for each thread. That object will have 
"thread scope" meaning it will be destructed only on thread termination. Note
that that means it can outlive destruction of the `ThreadLocal<T>` object.

To use it, simply create a `ThreadLocal<T>`, then call get() or upd() to get a 
readable or writable reference to the value of type `T` that is available for
the exclusive use of the current thread:
@code
    ThreadLocal<int> x;
    ...
    x.upd() = 5;
    assert(x.get() == 5);
@endcode

@warning You should avoid allocating `ThreadLocal` objects in single-threaded
code because the objects of type `T` have "thread scope"; they do not get 
destructed when the `ThreadLocal` object does. So in the single-threaded case
they will persist until program termination, creating a potential for leaks.
**/
template <class T>
class ThreadLocal {
public:
    /** Create a new `ThreadLocal<T>` object. This does not allocate any of
    the thread-local objects of type `T`; that is done from the individual
    threads when they first request such an object. **/
    ThreadLocal() {
        initialize();
    }

    /** Create a new `ThreadLocal<T>` object and provide a default value of type
    `T` to be used to initialize the thread-local objects when they are first
    allocated by the individual threads.

    @param defaultValue     The initial value that the objects of type `T` will
                            have when created by the individual threads.
    **/
    explicit ThreadLocal(const T& defaultValue) : m_defaultValue(defaultValue) {
        initialize();
    }

    /** Destructor deletes the thread local object but does not delete the
    individual thread-allocated objects of type T. Those are deleted only on
    thread termination. **/
    ~ThreadLocal() {
        pthread_key_delete(m_key);
    }

    /** Get a writable reference to the value of type `T` that was allocated
    for the current thread's exclusive use. **/
    T& upd() {
        T* value = reinterpret_cast<T*>(pthread_getspecific(m_key));
        if (!value) 
            value = createValue();
        return *value;
    }

    /** Get a const reference to the value of type `T` that was allocated
    for the current thread's exclusive use. **/
    const T& get() const {
        T* value = reinterpret_cast<T*>(pthread_getspecific(m_key));
        if (!value) 
            value = createValue();
        return *value;
    }

private:
    // This is a destructor for an object of type T that was allocated for this
    // thread's exclusive use. This will be called at thread termination, from
    // the same thread that allocated this object via the createValue() 
    // method below.
    static void cleanUpThreadLocalStorage(void* value) {
        T* t = reinterpret_cast<T*>(value);
        delete t;
    }

    void initialize() {
        pthread_key_create(&m_key, cleanUpThreadLocalStorage);
    }

    T* createValue() const {
        T* value = new T(m_defaultValue);
        pthread_setspecific(m_key, value);
        return value;
    }

    pthread_key_t   m_key;
    T               m_defaultValue;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_THREAD_LOCAL_H_
