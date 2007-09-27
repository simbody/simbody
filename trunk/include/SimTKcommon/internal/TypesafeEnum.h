#ifndef SimTK_SimTKCOMMON_TYPESAFE_ENUM_H_
#define SimTK_SimTKCOMMON_TYPESAFE_ENUM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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

#include "SimTKcommon.h"
#include <vector>

using std::vector;

namespace SimTK {

template <class T> class EnumSet;

/**
 * This class defines a typesafe enumeration.  It works like the "enum" keyword, but has several advantages:
 * 
 * - As the name implies, it is fully typesafe.  Where "enum" simply defines int values, the values defined by
 *   a TypesafeEnum are objects.
 * - It defines a static getAllValues() method which can be queried at runtime to programmatically determine
 *   the list of all possible values for the enumeration.  This is in contrast to "enum", where it is impossible
 *   to programmatically determine the list of values, or even how many allowed values there are.
 * - Because the values are objects, they can be extended in arbitrary ways.  You can define new methods or
 *   new metadata which should be associated with each value.
 * 
 * To create an enumeration, define a subclass of TypesafeEnum which is parameterized by itself.  For example:
 * 
 * <pre>
 * class Color : public TypesafeEnum<Color> {
 * public:
 *     static const Color Red;
 *     static const Color Green;
 *     static const Color Blue;
 * private:
 *     Color() : TypesafeEnum<Color>() {
 *     }
 * };
 * const Color Color::Red;
 * const Color Color::Green;
 * const Color Color::Blue;
 * </pre>
 * 
 * You can then define variables of the enumerated type and work with them exactly as you would expect:
 * 
 * <pre>
 * Color myColor = Color::Blue;
 * if (myColor == Color::Red) {...}
 * </pre>
 * 
 * And so on.  You can invoke getAllValues() to get the list of all possible values.  You also can invoke
 * getIndex() on any value of the enumerated type to find its index in the list:
 * 
 * <pre>
 * assert(Color::getAllValues()[Color::Red.getIndex()] == Color::Red);
 * </pre>
 * 
 * This class is designed to be used with the EnumSet class, which provides a convenient and efficient way
 * to represent a set of values of a particular enumerated type.
 */

template <class T>
class TypesafeEnum {
public:
    /**
     * Get the index of this value in the list returned by getAllValues().
     */
    int getIndex() const {
        return index;
    }
    /**
     * Get a list of all allowed values for this enumerated type.
     */
    static const vector<TypesafeEnum<T> >& getAllValues() {
        return updAllValues();
    }
    TypesafeEnum(const TypesafeEnum<T>& copy) : index(copy.index) {
    }
    T operator=(const TypesafeEnum<T>& copy) {
        index = copy.index;
    }
    bool operator==(const TypesafeEnum<T>& value) const {
        return (index == value.index);
    }
protected:
    TypesafeEnum() : index(updAllValues().size()) {
        SimTK_ASSERT_ALWAYS(index >= 0 && index <= 63, "Enum index is outside legal range.");
        int mask = 1<<index;
        SimTK_ASSERT_ALWAYS((allIndices&mask) == 0, "Multiple enumerated values have the same index.");
        allIndices |= mask;
        updAllValues().push_back(*this);
    }
private:
    int index;
    static int allIndices;
    void init() {
        SimTK_ASSERT_ALWAYS(index >= 0 && index <= 63, "Enum index is outside legal range.");
        int mask = 1<<index;
        SimTK_ASSERT_ALWAYS((allIndices&mask) == 0, "Multiple enumerated values have the same index.");
        allIndices |= mask;
        updAllValues().push_back(*this);
    }
    static vector<TypesafeEnum<T> >& updAllValues() {
        static vector<TypesafeEnum<T> > allValues;
        return allValues;
    }
};

template <class T>
int TypesafeEnum<T>::allIndices = 0;

/**
 * This class provides an efficient implementation of a set for storing values of an enumerated type
 * defined with TypesafeEnum.  The set is represented internally with bit flags, so storage, assignment,
 * and lookup are all extremely efficient.
 */

template <class T>
class EnumSet {
public:
    /**
     * Create an empty EnumSet.
     */
    EnumSet() : flags(0) {
    }
    /**
     * Create an EnumSet which contains a single value.
     */
    EnumSet(const TypesafeEnum<T>& value) : flags(mask(value)) {
    }
    /**
     * Create an EnumSet which contains the same values as another set.
     */
    EnumSet(const EnumSet<T>& set) : flags(set.flags) {
    }
    /**
     * Get the number of elements in this set.
     */
    int size() {
        // This method of counting the entries has linear performance in the number of entries.
        // Since enum sets will often contain only a small number of entries, this should be fast.
        
        int count = 0;
        int temp = flags;
        while (temp != 0) {
            ++count;
            temp &= temp-1;
        }
        return count;
    }
    /**
     * Determine whether this set contains a particular value.
     */
    bool contains(const TypesafeEnum<T>& value) {
        return ((flags & mask(value)) != 0);
    }
    /**
     * Determine whether this set contains all of the values in another set.
     */
    bool containsAll(const EnumSet<T>& set) {
        return ((flags & set.flags) == set.flags);
    }
    /**
     * Determine wheter this set contains any value which is in another set.
     */
    bool containsAny(const EnumSet<T>& set) {
        return ((flags & set.flags) != 0);
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator==(const EnumSet<T>& set) {
        return (flags == set.flags);
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator!=(const EnumSet<T>& set) {
        return (flags != set.flags);
    }
    EnumSet<T> operator=(const EnumSet<T>& set) {
        flags = set.flags;
        return *this;
    }
    EnumSet<T> operator+=(const TypesafeEnum<T>& value) {
        flags |= mask(value);
        return *this;
    }
    EnumSet<T> operator+=(const EnumSet<T>& set) {
        flags |= set.flags;
        return *this;
    }
    EnumSet<T> operator-=(const TypesafeEnum<T>& value) {
        flags &= ~mask(value);
        return *this;
    }
    EnumSet<T> operator-=(const EnumSet<T>& set) {
        flags &= ~set.flags;
        return *this;
    }
    EnumSet<T> operator+(const TypesafeEnum<T>& value) {
        EnumSet<T> temp(this);
        temp += value;
        return temp;
    }
    EnumSet<T> operator+(const EnumSet<T>& set) {
        EnumSet<T> temp(this);
        temp += set;
        return temp;
    }
    EnumSet<T> operator-(const TypesafeEnum<T>& value) {
        EnumSet<T> temp(this);
        temp -= value;
        return temp;
    }
    EnumSet<T> operator-(const EnumSet<T>& set) {
        EnumSet<T> temp(this);
        temp -= set;
        return temp;
    }
private:
    long mask(const TypesafeEnum<T>& value) const {
        return 1L<<value.getIndex();
    }
    long flags;
};

template <class T>
EnumSet<T> operator+(const TypesafeEnum<T>& value1, const TypesafeEnum<T>& value2) {
    EnumSet<T> temp(value1);
    temp += value2;
    return temp;
}

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_TYPESAFE_ENUM_H_
