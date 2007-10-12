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

using std::ostream;
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
 *     Color(int index, char* name) : TypesafeEnum<Color>(index, name) {
 *     }
 *     static void initValues() {
 *         new(&const_cast<Color&>(Red)) Color(0, "Red");
 *         new(&const_cast<Color&>(Green)) Color(1, "Green");
 *         new(&const_cast<Color&>(Blue)) Color(2, "Blue");
 *     }
 *     friend class TypesafeEnum<Color>;
 * };
 * 
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
        init();
        return index;
    }
    /**
     * Get the name of this value.
     */
    const char* getName() const {
        init();
        return name;
    }
    /**
     * Get a list of all allowed values for this enumerated type.
     */
    static const vector<TypesafeEnum<T> >& getAllValues() {
        init();
        return updAllValues();
    }
    TypesafeEnum(const TypesafeEnum<T>& copy) {
        init();
        index = copy.index;
        name = copy.name;
    }
    TypesafeEnum<T>& operator=(const TypesafeEnum<T>& copy) {
        init();
        index = copy.index;
        name = copy.name;
        return *this;
    }
    TypesafeEnum<T>* operator&() {
        init();
        return this;
    }
    bool operator==(const TypesafeEnum<T>& value) const {
        init();
        return (index == value.index);
    }
    bool operator!=(const TypesafeEnum<T>& value) const {
        init();
        return (index != value.index);
    }
    operator int() const {
        return getIndex();
    }
protected:
    TypesafeEnum() {
        init();
    }
    TypesafeEnum(int index, const char* name) : index(index), name(name) {
        SimTK_ASSERT_ALWAYS(index == updAllValues().size(), "Indices must be consecutive ints starting from 0.");
        int mask = 1<<index;
        updAllValues().push_back(*this);
    }
private:
    int index;
    const char* name;
    static vector<TypesafeEnum<T> >& updAllValues() {
        static vector<TypesafeEnum<T> > allValues;
        return allValues;
    }
    /**
     * Initialize the values of all enumerated constants.  This is invoked when a TypesafeEnum object is constructed,
     * or is referenced in any way.  This ensures that no one will ever see the static constants in an uninitialized
     * state, even if they are referenced by a static initializer in a different class that is initialized before
     * this one.
     */
    static void init() {
        static bool isInitialized = false;
        if (!isInitialized) {
            isInitialized = true;
            T::initValues();
        }
    }
};

template <class T>
ostream& operator<<(ostream& stream, const TypesafeEnum<T>& value) {
    return stream << value.getName();
}

/**
 * This class provides an efficient implementation of a set for storing values of an enumerated type
 * defined with TypesafeEnum.  The set is represented internally with bit flags, so storage, assignment,
 * and lookup are all extremely efficient.
 * 
 * This class supports all the standard bitwise operators, like &, |, ^, and ~.  This allows you to
 * manipulate sets exactly as if they were ints.  It also supports the + and - operators, which represent
 * the union and difference of two sets, respectively.
 * 
 * For example, if a method expects an EnumSet<Color> as an argument, you could pass any of the following
 * values:
 * 
 * <pre>
 * Color::Red
 * Color::Green+Color::Blue
 * EnumSet<Color>() // (an empty set)
 * ~EnumSet<Color>() // (the set of all possible values)
 * </pre>
 */

template <class T>
class EnumSet {
public:
    class iterator;
    /**
     * Create an empty EnumSet.
     */
    EnumSet() {
        init();
    }
    /**
     * Create an EnumSet which contains a single value.
     */
    EnumSet(const TypesafeEnum<T>& value) {
        init();
        flags[word(value)] = mask(value);
        numElements = 1;
    }
    /**
     * Create an EnumSet which contains the same values as another set.
     */
    EnumSet(const EnumSet<T>& set) {
        init();
        *this = set;
    }
    ~EnumSet() {
        delete[] flags;
    }
    /**
     * Get the number of elements in this set.
     */
    int size() const {
        if (numElements == -1) {
            numElements = 0;
            for (int i = 0; i < words; ++i) {
                int temp = flags[i];
                while (temp != 0) {
                    ++numElements;
                    temp &= temp-1;
                }
            }
        }
        return numElements;
    }
    /**
     * Determine whether this set contains a particular value.
     */
    bool contains(const TypesafeEnum<T>& value) const {
        return ((flags[word(value)] & mask(value)) != 0);
    }
    /**
     * Determine whether this set contains all of the values in another set.
     */
    bool containsAll(const EnumSet<T>& set) const {
        for (int i = 0; i < words; ++i)
            if ((flags[i] & set.flags[i]) != set.flags[i])
                return false;
        return true;
    }
    /**
     * Determine wheter this set contains any value which is in another set.
     */
    bool containsAny(const EnumSet<T>& set) const {
        for (int i = 0; i < words; ++i)
            if ((flags[i] & set.flags[i]) != 0)
                return true;
        return false;
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator==(const EnumSet<T>& set) const {
        for (int i = 0; i < words; ++i)
            if (flags[i] != set.flags[i])
                return false;
        return true;
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator!=(const EnumSet<T>& set) const {
        for (int i = 0; i < words; ++i)
            if (flags[i] != set.flags[i])
                return true;
        return false;
    }
    /**
     * Remove all elements from the set.
     */
    void clear() {
        for (int i = 0; i < words; ++i)
            flags[i] = 0;
        numElements = 0;
    }
    /**
     * Get an iterator pointing to the start of the set.
     */
    iterator begin() {
        iterator first(this, 0);
        return first;
    }
    /**
     * Get an iterator pointing to the end of the set.
     */
    iterator end() {
        iterator last(this, TypesafeEnum<T>::getAllValues().size());
        return last;
    }
    EnumSet<T>& operator=(const EnumSet<T>& set) {
        for (int i = 0; i < words; ++i)
            flags[i] = set.flags[i];
        numElements = set.numElements;
        return *this;
    }
    EnumSet<T>& operator+=(const TypesafeEnum<T>& value) {
        flags[word(value)] |= mask(value);
        numElements = -1;
        return *this;
    }
    EnumSet<T>& operator+=(const EnumSet<T>& set) {
        for (int i = 0; i < words; ++i)
            flags[i] |= set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumSet<T>& operator-=(const TypesafeEnum<T>& value) {
        flags[word(value)] &= ~mask(value);
        numElements = -1;
        return *this;
    }
    EnumSet<T>& operator-=(const EnumSet<T>& set) {
        for (int i = 0; i < words; ++i)
            flags[i] &= ~set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumSet<T> operator+(const TypesafeEnum<T>& value) const {
        EnumSet<T> temp(this);
        temp += value;
        return temp;
    }
    EnumSet<T> operator+(const EnumSet<T>& set) const {
        EnumSet<T> temp(*this);
        temp += set;
        return temp;
    }
    EnumSet<T> operator-(const TypesafeEnum<T>& value) const {
        EnumSet<T> temp(this);
        temp -= value;
        return temp;
    }
    EnumSet<T> operator-(const EnumSet<T>& set) const {
        EnumSet<T> temp(*this);
        temp -= set;
        return temp;
    }
    EnumSet<T>& operator|=(const EnumSet<T>& set) {
        for (int i = 0; i < words; ++i)
            flags[i] |= set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumSet<T> operator|(const EnumSet<T>& set) const {
        EnumSet<T> temp(*this);
        temp += set;
        return temp;
    }
    EnumSet<T>& operator&=(const EnumSet<T>& set) {
        for (int i = 0; i < words; ++i)
            flags[i] &= set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumSet<T> operator&(const EnumSet<T>& set) const {
        EnumSet<T> temp(*this);
        temp &= set;
        return temp;
    }
    EnumSet<T>& operator^=(const EnumSet<T>& set) {
        for (int i = 0; i < words; ++i)
            flags[i] ^= set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumSet<T> operator^(const EnumSet<T>& set) const {
        EnumSet<T> temp(*this);
        temp ^= set;
        return temp;
    }
    EnumSet<T> operator~() const {
        EnumSet<T> temp(*this);
        for (int i = 0; i < words-1; ++i)
            temp.flags[i] = -1-temp.flags[i];
        temp.flags[words-1] = (1<<T::getAllValues().size()%BITS_PER_WORD)-1-flags[words-1];
        temp.numElements = -1;
        return temp;
    }
private:
    /**
     * Initialize the set.
     */
    void init() {
        numElements = -1;
        words = (T::getAllValues().size()+BITS_PER_WORD-1)/BITS_PER_WORD;
        flags = new int[words];
        for (int i = 0; i < words; ++i)
            flags[i] = 0;
    }
    /**
     * Given an element of the enumeration, return a mask in which the bit for that element is set.
     * The mask is for the word that contains the specified element.
     */
    int mask(const TypesafeEnum<T>& value) const {
        return 1<<(value.getIndex()%BITS_PER_WORD);
    }
    /**
     * Given an element of the enumeration, return the index of the word that contains that element.
     */
    int word(const TypesafeEnum<T>& value) const {
        return value.getIndex()/BITS_PER_WORD;
    }
    int* flags;
    short words;
    mutable short numElements;
    static const int BITS_PER_WORD = 32;
};

template <class T>
EnumSet<T> operator+(const TypesafeEnum<T>& value1, const TypesafeEnum<T>& value2) {
    EnumSet<T> temp(value1);
    temp += value2;
    return temp;
}

/**
 * This class provides an interface for iterating over the content of an EnumSet.
 */

template <class T>
class EnumSet<T>::iterator {
public:
    TypesafeEnum<T> operator*() {
        assert (index >= 0 && index < TypesafeEnum<T>::getAllValues().size());
        return TypesafeEnum<T>::getAllValues()[index];
    }
    iterator operator++() {
        assert (index < TypesafeEnum<T>::getAllValues().size());
        ++index;
        findNextElement();
        return *this;
    }
    iterator operator++(int) {
        assert (index < TypesafeEnum<T>::getAllValues().size());
        iterator current = *this;
        ++index;
        findNextElement();
        return current;
    }
    bool operator==(const iterator& iter) {
        return (set == iter.set && index == iter.index);
    }
    bool operator!=(const iterator& iter) {
        return (set != iter.set || index != iter.index);
    }
private:
    iterator(EnumSet<T>* set, int index) : set(set), index(index) {
        findNextElement();
    }
    void findNextElement() {
        while (index < TypesafeEnum<T>::getAllValues().size() && !set->contains(TypesafeEnum<T>::getAllValues()[index]))
            ++index;
    }
    int index;
    EnumSet<T>* set;
    friend class EnumSet<T>;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_TYPESAFE_ENUM_H_
