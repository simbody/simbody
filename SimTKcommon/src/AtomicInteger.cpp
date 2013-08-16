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

#include "SimTKcommon/internal/AtomicInteger.h"
#include "gmx_atomic.h"

using namespace SimTK;

AtomicInteger::AtomicInteger() {
    atomic = new gmx_atomic_t();
    gmx_atomic_set(reinterpret_cast<gmx_atomic_t*>(atomic), 0);
}

AtomicInteger::AtomicInteger(int value) {
    atomic = new gmx_atomic_t();
    gmx_atomic_set(reinterpret_cast<gmx_atomic_t*>(atomic), value);
}

AtomicInteger::~AtomicInteger() {
    delete reinterpret_cast<gmx_atomic_t*>(atomic);
}

AtomicInteger& AtomicInteger::operator=(int value) {
    gmx_atomic_set(reinterpret_cast<gmx_atomic_t*>(atomic), value);
    return *this;
}

AtomicInteger::operator int() const {
    return gmx_atomic_read(reinterpret_cast<gmx_atomic_t*>(atomic));
}

int AtomicInteger::operator++() {
    return gmx_atomic_add_return(reinterpret_cast<gmx_atomic_t*>(atomic), 1);
}

int AtomicInteger::operator++(int) {
    return gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomic), 1);
}

int AtomicInteger::operator--() {
    return gmx_atomic_add_return(reinterpret_cast<gmx_atomic_t*>(atomic), -1);
}

int AtomicInteger::operator--(int) {
    return gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomic), -1);
}

AtomicInteger& AtomicInteger::operator+=(int value) {
    gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomic), value);
    return *this;
}

AtomicInteger& AtomicInteger::operator-=(int value) {
    gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomic), -value);
    return *this;
}

AtomicInteger& AtomicInteger::operator*=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue*value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator/=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue/value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator%=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue%value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator&=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue&value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator|=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue|value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator^=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue^value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator<<=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue<<value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

AtomicInteger& AtomicInteger::operator>>=(int value) {
    gmx_atomic_t* a = reinterpret_cast<gmx_atomic_t*>(atomic);
    int oldvalue, newvalue;
    do {
        oldvalue = gmx_atomic_read(a);
        newvalue = oldvalue>>value;
    } while (gmx_atomic_cmpxchg(a, oldvalue, newvalue) != oldvalue);
    return *this;
}

bool AtomicInteger::operator==(int value) const {
    return (gmx_atomic_read(reinterpret_cast<gmx_atomic_t*>(atomic)) == value);
}

bool AtomicInteger::operator!=(int value) const {
    return (gmx_atomic_read(reinterpret_cast<gmx_atomic_t*>(atomic)) != value);
}

std::ostream& SimTK::operator<<(std::ostream& stream, const AtomicInteger& value) {
    int i = value;
    return stream << i;
}
