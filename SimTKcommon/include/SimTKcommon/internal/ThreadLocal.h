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

#include <pthread.h>
#include <map>
#include <set>

namespace SimTK {

static std::map<void*, pthread_key_t> instanceMap;
static std::map<pthread_key_t, std::set<void*> > keyInstances;
static pthread_mutex_t keyLock = PTHREAD_MUTEX_INITIALIZER;

template <class T>
static void cleanUpThreadLocalStorage(void* value) {

    // Delete the value.

    T* t = reinterpret_cast<T*>(value);
    delete t;

    // Remove it from the set of values needing to be deleted.

    pthread_mutex_lock(&keyLock);
    pthread_key_t key = instanceMap[value];
    instanceMap.erase(value);
    if (keyInstances.find(key) != keyInstances.end())
        keyInstances.find(key)->second.erase(value);
    pthread_mutex_unlock(&keyLock);

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
        this->initialize();
    }
    /**
     * Create a new ThreadLocal variable.
     *
     * @param defaultValue the initial value which the variable will have on each thread
     */
    ThreadLocal(const T& defaultValue) : defaultValue(defaultValue) {
        this->initialize();
    }
    ~ThreadLocal() {

        // Delete the key.

        pthread_key_delete(key);

        // Once the key is deleted, cleanUpThreadLocalStorage() will no longer be called, so delete
        // all instances now.

        pthread_mutex_lock(&keyLock);
        std::set<void*>& instances = keyInstances[key];
        for (std::set<void*>::const_iterator iter = instances.begin(); iter != instances.end(); ++iter) {
            instanceMap.erase(*iter);
            delete (T*) *iter;
        }
        keyInstances.erase(key);
        pthread_mutex_unlock(&keyLock);
    }
    /**
     * Get a reference to the value for the current thread.
     */
    T& upd() {
        T* value = reinterpret_cast<T*>(pthread_getspecific(key));
        if (value == NULL)
            return createValue();
        return *value;
    }
    /**
     * Get a const reference to the value for the current thread.
     */
    const T& get() const {
        T* value = reinterpret_cast<T*>(pthread_getspecific(key));
        if (value == NULL)
            return createValue();
        return *value;
    }
private:
    void initialize() {
        pthread_key_create(&key, cleanUpThreadLocalStorage<T>);
        pthread_mutex_lock(&keyLock);
        keyInstances[key] = std::set<void*>();
        pthread_mutex_unlock(&keyLock);
    }
    T& createValue() const {
        T* value = new T(defaultValue);
        pthread_setspecific(key, value);
        pthread_mutex_lock(&keyLock);
        instanceMap[value] = key;
        keyInstances[key].insert(value);
        pthread_mutex_unlock(&keyLock);
        return *value;
    }
    pthread_key_t key;
    T defaultValue;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_THREAD_LOCAL_H_
