#ifndef SimTK_SimTKCOMMON_ATOMIC_INTEGER_H_
#define SimTK_SimTKCOMMON_ATOMIC_INTEGER_H_

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

#include "SimTKcommon/internal/common.h"

namespace SimTK {

/**
 * This class functions exactly like an int, except that the following operators are atomic:
 * ++, --, +=, -=, *=, /=, %=, &=, |=, ^=, <<=, and >>=.  For example, suppose myInt is an AtomicInteger
 * that initially has the value 5.  If two threads both evaluate the expression ++myInt, it is guaranteed
 * that one thread will get the value 6 and the other will get the value 7, and myInt will
 * equal 7 afterward.  This would not be true for an ordinary int.
 * 
 * On most processors, this form of thread-safety can be implemented in a lightweight way
 * which is much faster than acquiring a lock.  When possible, this class uses these mechanisms
 * to achieve maximum efficiency.  On platforms that do not support atomic operations directly
 * it uses locking, which is slower but still guaranteed to produce a correct result.
 */

class SimTK_SimTKCOMMON_EXPORT AtomicInteger {
public:
    AtomicInteger();
    AtomicInteger(int value);
    ~AtomicInteger();
    AtomicInteger& operator=(int value);
    operator int() const;
    int operator++();
    int operator++(int);
    int operator--();
    int operator--(int);
    AtomicInteger& operator+=(int value);
    AtomicInteger& operator-=(int value);
    AtomicInteger& operator*=(int value);
    AtomicInteger& operator/=(int value);
    AtomicInteger& operator%=(int value);
    AtomicInteger& operator&=(int value);
    AtomicInteger& operator|=(int value);
    AtomicInteger& operator^=(int value);
    AtomicInteger& operator<<=(int value);
    AtomicInteger& operator>>=(int value);
    bool operator==(int value) const;
    bool operator!=(int value) const;
private:
    void* atomic;
};

std::ostream& operator<<(std::ostream& stream, const AtomicInteger& value);

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_ATOMIC_INTEGER_H_
