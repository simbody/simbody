#ifndef SimTK_SimTKCOMMON_ATOMIC_INTEGER_H_
#define SimTK_SimTKCOMMON_ATOMIC_INTEGER_H_

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

#include "common.h"

namespace SimTK {

/**
 * This class functions exactly like an int, except that the following operators are atomic:
 * ++, --, +=, -=, *=, and /=.  For example, suppose myInt is an AtomicInteger that initially
 * has the value 5.  If two threads both evaluate the expression ++myInt, it is guaranteed
 * that one thread will get the value 6 and the other will get the value 7, and myInt will
 * equal 7 afterward.  This would not be true for an ordinary int.
 * 
 * On most processors, this form of thread-safety can be implemented in a lightweight way
 * which is much faster than acquiring a lock.  When possible, this class uses these mechanisms
 * to achieve maximum efficiency.  On platforms that do not support atomic operations directly,
 * it uses locking which is slower, but still guaranteed to produce a correct result.
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
    bool operator==(int value) const;
    bool operator!=(int value) const;
private:
    void* atomic;
};

std::ostream& operator<<(std::ostream& stream, const AtomicInteger& value);

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_ATOMIC_INTEGER_H_
