#ifndef SimTK_SimTKCOMMON_THREAD_LOCAL_H_
#define SimTK_SimTKCOMMON_THREAD_LOCAL_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include <pthread.h>

namespace SimTK {

template <class T>
static void cleanUpThreadLocalStorage(void* value) {
    if (value != NULL) {
        T* t = reinterpret_cast<T*>(value);
        delete t;
    }
}

/**
 * This class represents a "thread local" variable: one which has a different value on each thread.
 * This is useful in many situations when writing multithreaded code.  For example, it can be used
 * as temporary workspace for calculations.  If a single workspace object were created, all access
 * to it would need to be synchronized to prevent threads from overwriting each other's values.
 * Using a ThreadLocal instead means that a separate workspace object will automatically be created
 * for each thread.
 * 
 * To use it, simply create a ThreadLocal, then call get() or upd() to get a readable or writable
 * reference to the value for the current thread:
 * 
 * <pre>
 * ThreadLocal<int> x;
 * ...
 * x.upd() = 5;
 * assert(x.get() == 5);
 * </pre>
 */

template <class T>
class ThreadLocal {
public:
    /**
     * Create a new ThreadLocal variable.
     */
    ThreadLocal() {
        pthread_key_create(&key, cleanUpThreadLocalStorage<T>);
    }
    /**
     * Create a new ThreadLocal variable.
     * 
     * @param defaultValue the initial value which the variable will have on each thread
     */
    ThreadLocal(const T& defaultValue) : defaultValue(defaultValue) {
        pthread_key_create(&key, cleanUpThreadLocalStorage<T>);
    }
    ~ThreadLocal() {
        pthread_key_delete(key);
    }
    /**
     * Get a reference to the value for the current thread.
     */
    T& upd() {
        T* value = reinterpret_cast<T*>(pthread_getspecific(key));
        if (value == NULL) {
            value = new T(defaultValue);
            pthread_setspecific(key, value);
        }
        return *value;
    }
    /**
     * Get a const reference to the value for the current thread.
     */
    const T& get() const {
        T* value = reinterpret_cast<T*>(pthread_getspecific(key));
        if (value == NULL) {
            value = new T(defaultValue);
            pthread_setspecific(key, value);
        }
        return *value;
    }
private:
    pthread_key_t key;
    T defaultValue;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_THREAD_LOCAL_H_
