#ifndef SimTK_SimTKCOMMON_ENUMERATION_H_
#define SimTK_SimTKCOMMON_ENUMERATION_H_

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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include <vector>

namespace SimTK {

template <class T> class EnumerationSet;

/**
 * This class defines an enumerated type.  It works like the "enum" keyword, but has several advantages:
 * 
 * - It is fully typesafe.  If an Enumeration defines a list of allowed values, it is absolutely guaranteed that any
 *   variable of that type will contain one of the allowed values, and no other type can be assigned to it.  In contrast,
 *   "enum" simply defines integer values, and it is possible for a variable of that type to contain any integer, not
 *   just one of the allowed values.
 * - It defines static methods which can be queried at runtime to programmatically determine the list of possible values
 *   for the enumeration.  This is in contrast to "enum", where it is impossible to programmatically determine the list
 *   of values, or even how many allowed values there are.
 * - Because the values are objects, they can be extended in arbitrary ways.  You can define new methods or
 *   new metadata which should be associated with each value.
 * 
 * One drawback is that there is no way to provide default construction for objects of Enumeration type, because
 * default construction is necessary internally to get the static values initialized properly. Consequently
 * you must declare the default constructor as private in your derived class.
 *
 * To create an enumeration, define a subclass of Enumeration which is parameterized by itself.  For example:
 * 
 * <pre>
 * class Color : public Enumeration<Color> {
 * public:
 *     enum Index {RedIndex = 0, GreenIndex = 1, BlueIndex = 2};
 *     static const Color Red;
 *     static const Color Green;
 *     static const Color Blue;
 * private:
 *     Color();
 *     Color(const Color& thisElement, int index, char* name);
 *     static void initValues();
 *     friend class Enumeration<Color>;
 * };
 * 
 * const Color Color::Red;
 * const Color Color::Green;
 * const Color Color::Blue;
 * 
 * Color::Color() : Enumeration<Color>() {
 * }
 * 
 * Color::Color(const Color& thisElement, int index, char* name) : Enumeration<Color>(thisElement, index, name) {
 * }
 * 
 * void Color::initValues() {
 *     new(&const_cast<Color&>(Red)) Color(Red, RedIndex, "Red");
 *     new(&const_cast<Color&>(Green)) Color(Green, GreenIndex, "Green");
 *     new(&const_cast<Color&>(Blue)) Color(Blue, BlueIndex, "Blue");
 * }
 * </pre>
 * 
 * You can then define variables of the enumerated type and work with them exactly as you would expect:
 * 
 * <pre>
 * Color myColor = Color::Blue;
 * if (myColor == Color::Red) {...}
 * </pre>
 * 
 * Enumerations can be used in switch statements, but the cases must be value indices rather than the values
 * themselves:
 * 
 * <pre>
 * switch (myColor) {
 *   case Color::RedIndex:
 *     ...
 * }
 * </pre>
 * 
 * You can invoke size() to determine the number of allowed values and getValue() to get the value with a
 * particular index.  Alternatively, you can iterate over the allowed values with begin() and end().  You
 * also can invoke getIndex() on any value of the enumerated type to find its index in the list:
 * 
 * <pre>
 * assert(Color::getValue(Color::Red.getIndex()) == Color::Red);
 * </pre>
 * 
 * This class is designed to be used with the EnumerationSet class, which provides a convenient and efficient way
 * to represent a set of values of a particular enumerated type.
 */

template <class T>
class Enumeration {
public:
    class iterator;
    /**
     * Get the index of this value in the list returned by getValue().
     */
    int getIndex() const {
        init();
        return index;
    }
    /**
     * Get the name of this value.
     */
    const std::string getName() const {
        init();
        return name;
    }
    /**
     * Get the total number of allowed values for this enumerated type.
     */
    static int size() {
        return (int)updAllValues().size();
    }
    /**
     * Get the enumerated value with a particular index.
     */
    static const T& getValue(int index) {
        return *updAllValues()[index];
    }
    /**
     * Get an iterator pointing to the start of the set of all possible values.
     */
    static iterator begin() {
        iterator first(0);
        return first;
    }
    /**
     * Get an iterator pointing to the end of the set of all possible values.
     */
    static iterator end() {
        iterator last(size());
        return last;
    }
    Enumeration(const Enumeration<T>& copy) {
        init();
        index = copy.index;
        name = copy.name;
    }
    Enumeration<T>& operator=(const Enumeration<T>& copy) {
        init();
        index = copy.index;
        name = copy.name;
        return *this;
    }
    // Address operator should usually return a Stage*, not a Enumeration<Stage>*
    // The previous Enumeration<Stage>* return value was causing a problem in 
    // boost.python headers, in a statment equivalent to:
    //    "Stage* const p = &const_cast<Stage&>((Stage const&)stage);"
    // Which looks innocent, except that it fails to compile when that first "&" 
    // generated an Enumeration<Stage>* instead of a Stage*.  Wrapping now works with
    // this method returning a T*.
    T* operator&() {
        init();
        return static_cast<T*>(this);
    }
    bool operator==(const Enumeration<T>& value) const {
        init();
        return (index == value.index);
    }
    bool operator!=(const Enumeration<T>& value) const {
        init();
        return (index != value.index);
    }
    operator int() const {
        return getIndex();
    }
    EnumerationSet<T> operator|(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp |= set;
        return temp;
    }
    EnumerationSet<T> operator&(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp &= set;
        return temp;
    }
    EnumerationSet<T> operator^(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp ^= set;
        return temp;
    }
    /**
     * Get the next enumerated value in order of their indices.
     */
    Enumeration<T> operator++() {
        assert (index < size()-1);
        ++index;
        name = getValue(index).name;
        return *this;
    }
    /**
     * Get the next enumerated value in order of their indices.
     */
    Enumeration<T> operator++(int) {
        assert (index < size()-1);
        Enumeration<T> current = *this;
        ++index;
        name = getValue(index).name;
        return current;
    }
    /**
     * Get the previous enumerated value in order of their indices.
     */
    Enumeration<T> operator--() {
        assert (index > 0);
        --index;
        name = getValue(index).name;
        return *this;
    }
    /**
     * Get the previous enumerated value in order of their indices.
     */
    Enumeration<T> operator--(int) {
        assert (index > 0);
        Enumeration<T> current = *this;
        --index;
        name = getValue(index).name;
        return current;
    }
protected:
    // Default construction must not do anything to memory, because the static members
    // may already have been initialized before they get default constructed. If you
    // touch any data here, say by initializing "index" to -1, you may be wiping out a
    // static member constants. The constant will not then get fixed because the init() 
    // call will think it is already done.
    Enumeration() {
        init();
    }
    Enumeration(const T& thisElement, int index, const char* name) : index(index), name(name) {
        SimTK_ASSERT_ALWAYS(index == updAllValues().size(), "Indices must be consecutive ints starting from 0.");
        updAllValues().push_back(&thisElement);
    }
    static std::vector<const T*>& updAllValues() {
        static std::vector<const T*> allValues;
        return allValues;
    }
private:
    int index;
    const char* name;
    /**
     * Initialize the values of all enumerated constants.  This is invoked when an Enumeration object is constructed,
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
std::ostream& operator<<(std::ostream& stream, const Enumeration<T>& value) {
    return stream << value.getName();
}

/**
 * This class provides an interface for iterating over the set of all possible enumerated values.
 */

template <class T>
class Enumeration<T>::iterator {
public:
    Enumeration<T> operator*() {
        assert (index >= 0 && index < Enumeration<T>::size());
        return Enumeration<T>::getValue(index);
    }
    iterator operator++() {
        assert (index < Enumeration<T>::size());
        ++index;
        return *this;
    }
    iterator operator++(int) {
        assert (index < Enumeration<T>::size());
        iterator current = *this;
        ++index;
        return current;
    }
    bool operator==(const iterator& iter) {
        return (index == iter.index);
    }
    bool operator!=(const iterator& iter) {
        return (index != iter.index);
    }
private:
    iterator(int index) : index(index) {
    }
    int index;
    friend class Enumeration<T>;
};

/**
 * This class provides an efficient implementation of a set for storing values of an enumerated type
 * defined with Enumeration.  The set is represented internally with bit flags, so storage, assignment,
 * and lookup are all extremely efficient.
 * 
 * This class supports all the standard bitwise operators, like &, |, ^, and ~.  This allows you to
 * manipulate sets exactly as if they were ints.  It also supports the - operator, which represents
 * the difference between two sets.
 * 
 * For example, if a method expects an EnumerationSet<Color> as an argument, you could pass any of the following
 * values:
 * 
 * <pre>
 * Color::Red                   // a set containing Red
 * Color::Green | Color::Blue   // a set containing Green and Blue
 * EnumerationSet<Color>()      // an empty set
 * ~EnumerationSet<Color>()     // the set of all possible values
 * </pre>
 */

template <class T>
class EnumerationSet {
public:
    class iterator;
    /**
     * Create an empty EnumerationSet.
     */
    EnumerationSet() {
        rep = new EnumerationSetRep();
    }
    /**
     * Create an EnumerationSet which contains a single value.
     */
    EnumerationSet(const Enumeration<T>& value) {
        rep = new EnumerationSetRep(value);
    }
    /**
     * Create an EnumerationSet which contains the same values as another set.
     */
    EnumerationSet(const EnumerationSet<T>& set) {
        rep = new EnumerationSetRep(*set.rep);
    }
    ~EnumerationSet() {
        delete rep;
    }
    /**
     * Get the number of elements in this set.
     */
    int size() const {
        return rep->size();
    }
    /**
     * Check whether this set is empty.
     */
    bool empty() const {
        return (rep->size() == 0);
    }
    /**
     * Determine whether this set contains a particular value.
     */
    bool contains(const Enumeration<T>& value) const {
        return rep->contains(value);
    }
    /**
     * Determine whether this set contains all of the values in another set.
     */
    bool containsAll(const EnumerationSet<T>& set) const {
        return rep->containsAll(*set.rep);
    }
    /**
     * Determine wheter this set contains any value which is in another set.
     */
    bool containsAny(const EnumerationSet<T>& set) const {
        return rep->containsAny(*set.rep);
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator==(const EnumerationSet<T>& set) const {
        return *rep == *set.rep;
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator!=(const EnumerationSet<T>& set) const {
        return *rep != *set.rep;
    }
    /**
     * Remove all elements from the set.
     */
    void clear() {
        rep->clear();
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
        iterator last(this, Enumeration<T>::size());
        return last;
    }
    EnumerationSet<T>& operator=(const EnumerationSet<T>& set) {
        *rep = *set.rep;
        return *this;
    }
    EnumerationSet<T>& operator-=(const Enumeration<T>& value) {
        *rep -= value;
        return *this;
    }
    EnumerationSet<T>& operator-=(const EnumerationSet<T>& set) {
        *rep -= *set.rep;
        return *this;
    }
    EnumerationSet<T> operator-(const Enumeration<T>& value) const {
        EnumerationSet<T> temp(*this);
        temp -= value;
        return temp;
    }
    EnumerationSet<T> operator-(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp -= set;
        return temp;
    }
    EnumerationSet<T>& operator|=(const EnumerationSet<T>& set) {
        *rep |= *set.rep;
        return *this;
    }
    EnumerationSet<T> operator|(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp |= set;
        return temp;
    }
    EnumerationSet<T>& operator&=(const EnumerationSet<T>& set) {
        *rep &= *set.rep;
        return *this;
    }
    EnumerationSet<T> operator&(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp &= set;
        return temp;
    }
    EnumerationSet<T>& operator^=(const EnumerationSet<T>& set) {
        *rep ^= *set.rep;
        return *this;
    }
    EnumerationSet<T> operator^(const EnumerationSet<T>& set) const {
        EnumerationSet<T> temp(*this);
        temp ^= set;
        return temp;
    }
    EnumerationSet<T> operator~() const {
        EnumerationSet<T> temp(*this);
        temp.rep->invert();
        return temp;
    }
private:
    class EnumerationSetRep;
    EnumerationSetRep* rep;
};

/**
 * This class is the internal implementation of EnumerationSet.
 */

template <class T>
class EnumerationSet<T>::EnumerationSetRep {
public:
    /**
     * Create an empty EnumerationSet.
     */
    EnumerationSetRep() {
        init();
    }
    /**
     * Create an EnumerationSet which contains a single value.
     */
    EnumerationSetRep(const Enumeration<T>& value) {
        init();
        flags[word(value)] = mask(value);
        numElements = 1;
    }
    /**
     * Create an EnumerationSet which contains the same values as another set.
     */
    EnumerationSetRep(const EnumerationSetRep& set) {
        init();
        *this = set;
    }
    ~EnumerationSetRep() {
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
    bool contains(const Enumeration<T>& value) const {
        return ((flags[word(value)] & mask(value)) != 0);
    }
    /**
     * Determine whether this set contains all of the values in another set.
     */
    bool containsAll(const EnumerationSetRep& set) const {
        for (int i = 0; i < words; ++i)
            if ((flags[i] & set.flags[i]) != set.flags[i])
                return false;
        return true;
    }
    /**
     * Determine wheter this set contains any value which is in another set.
     */
    bool containsAny(const EnumerationSetRep& set) const {
        for (int i = 0; i < words; ++i)
            if ((flags[i] & set.flags[i]) != 0)
                return true;
        return false;
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator==(const EnumerationSetRep& set) const {
        for (int i = 0; i < words; ++i)
            if (flags[i] != set.flags[i])
                return false;
        return true;
    }
    /**
     * Determine whether this set has identical contents to another one.
     */
    bool operator!=(const EnumerationSetRep& set) const {
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
    EnumerationSetRep& operator=(const EnumerationSetRep& set) {
        for (int i = 0; i < words; ++i)
            flags[i] = set.flags[i];
        numElements = set.numElements;
        return *this;
    }
    EnumerationSetRep& operator-=(const Enumeration<T>& value) {
        flags[word(value)] &= ~mask(value);
        numElements = -1;
        return *this;
    }
    EnumerationSetRep& operator-=(const EnumerationSetRep& set) {
        for (int i = 0; i < words; ++i)
            flags[i] &= ~set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumerationSetRep& operator|=(const EnumerationSetRep& set) {
        for (int i = 0; i < words; ++i)
            flags[i] |= set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumerationSetRep& operator&=(const EnumerationSetRep& set) {
        for (int i = 0; i < words; ++i)
            flags[i] &= set.flags[i];
        numElements = -1;
        return *this;
    }
    EnumerationSetRep& operator^=(const EnumerationSetRep& set) {
        for (int i = 0; i < words; ++i)
            flags[i] ^= set.flags[i];
        numElements = -1;
        return *this;
    }
    void invert() const {
        for (int i = 0; i < words-1; ++i)
            flags[i] = -1-flags[i];
        flags[words-1] = (1<<T::size()%BITS_PER_WORD)-1-flags[words-1];
        numElements = -1;
    }
private:
    /**
     * Initialize the set.
     */
    void init() {
        numElements = -1;
        words = (T::size()+BITS_PER_WORD-1)/BITS_PER_WORD;
        flags = new int[words];
        for (int i = 0; i < words; ++i)
            flags[i] = 0;
    }
    /**
     * Given an element of the enumeration, return a mask in which the bit for that element is set.
     * The mask is for the word that contains the specified element.
     */
    int mask(const Enumeration<T>& value) const {
        return 1<<(value.getIndex()%BITS_PER_WORD);
    }
    /**
     * Given an element of the enumeration, return the index of the word that contains that element.
     */
    int word(const Enumeration<T>& value) const {
        return value.getIndex()/BITS_PER_WORD;
    }
    int* flags;
    short words;
    mutable short numElements;
    static const int BITS_PER_WORD = 8*sizeof(int);
};

/**
 * This class provides an interface for iterating over the content of an EnumerationSet.
 */

template <class T>
class EnumerationSet<T>::iterator {
public:
    Enumeration<T> operator*() {
        assert (index >= 0 && index < Enumeration<T>::size());
        return Enumeration<T>::getValue(index);
    }
    iterator operator++() {
        assert (index < Enumeration<T>::size());
        ++index;
        findNextElement();
        return *this;
    }
    iterator operator++(int) {
        assert (index < Enumeration<T>::size());
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
    iterator(EnumerationSet<T>* set, int index) : set(set), index(index) {
        findNextElement();
    }
    void findNextElement() {
        while (index < Enumeration<T>::size() && !set->contains(Enumeration<T>::getValue(index)))
            ++index;
    }
    int index;
    EnumerationSet<T>* set;
    friend class EnumerationSet<T>;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_ENUMERATION_H_
