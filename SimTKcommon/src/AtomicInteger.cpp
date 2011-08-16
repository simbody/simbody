/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
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
