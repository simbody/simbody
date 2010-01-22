/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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

/*
 * These are regression tests for the SimTK::Array_<T,X> class.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include <vector>
#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

SimTK_DEFINE_UNIQUE_INDEX_TYPE(TestIx);

class SmallIx {
public:
    SmallIx() : ix(0xff) {}
    explicit SmallIx(unsigned char i) : ix(i) {}

    SmallIx& operator++() 
    {   assert(ix<max_size); ++ix; return *this;}
    SmallIx operator++(int) 
    {   assert(ix<max_size); const SmallIx x=*this; ++ix; return x;}
    SmallIx& operator--() 
    {   assert(ix>0); --ix; return *this;}
    SmallIx operator--(int) 
    {   assert(ix>0); const SmallIx x=*this; ++ix; return x;}
    operator unsigned char() const 
    {   return ix;}

    typedef unsigned char index_type;
    typedef unsigned char size_type;
    typedef signed char   difference_type;
    static const size_type max_size = 4;
    static const char* index_name() {return "SmallIx";}
private:
    unsigned char ix;
};

class Counter {
public:
    Counter() : count(0) {}
    Counter& operator=(int i) {count=i; return *this;}
    Counter& operator++() {++count; return *this;}
    Counter operator++(int) {const Counter c=*this; ++count; return c;}
    Counter& reset() {count=0; return *this;}
    operator int() const {return count;}
private:
    mutable int count;
};
inline std::ostream&
operator<<(std::ostream& o, const Counter& c) {
    return o << (int)c;
} 

template <class T>
struct Count {
    Count() {++defCtor;}
    Count(const T& t) : val(t) {++initCtor;}
    Count(const Count& c) : val(c.val) {++copyCtor;}
    ~Count() {++dtor;}

    static void dumpCounts(const char* msg) {
        cout << msg << ":";
        cout << " def=" << Count<int>::defCtor;
        cout << " init=" << Count<int>::initCtor;
        cout << " copy=" << Count<int>::copyCtor;
        cout << " dtor=" << Count<int>::dtor;
        cout << endl;
    }

    T val;

    static void reset() {defCtor=initCtor=copyCtor=dtor=0;}
    static Counter defCtor;
    static Counter initCtor;
    static Counter copyCtor;
    static Counter dtor;
};
template <class T> inline std::ostream&
operator<<(std::ostream& o, const Count<T>& c) {
    return o << c.val;
} 
template <class T> Counter Count<T>::defCtor;
template <class T> Counter Count<T>::initCtor;
template <class T> Counter Count<T>::copyCtor;
template <class T> Counter Count<T>::dtor;

void testConstruction() {
    const int data[] = {5,3,-2,27,9};
    const char uchar[] = {'f','i','t','z'};
    Array_<int> nothing;
    Array_<int> def(5);
    Array_<int> intWithInt(data, data+5);
    Array_<char> charWithChar(uchar, uchar+4);
    Array_<int>  intWithChar(uchar, uchar+4);
    Array_<char>  charWithInt(data, data+5);
    cout << "nothing=" << nothing << endl;
    cout << "def=" << def << endl;
    cout << "intWithInt=" << intWithInt << endl;
    cout << "charWithChar=" << charWithChar << endl;
    cout << "intWithChar=" << intWithChar << endl;
    cout << "charWithInt=" << charWithInt << endl;

    Array_< Count<int> > cint(data, data+5);
    Count<int>::dumpCounts("cint(data,data+5)");
    Count<int>::reset();

    const Count<int> counts[] = {3,4,5};
    Count<int>::reset();
    Array_< Count<int> > ccnt(counts, counts+3);
    Count<int>::dumpCounts("ccnt(counts,counts+3)");
    Count<int>::reset();

    Array_< Count<int> > cint2(cint);
    Count<int>::dumpCounts("cint2(cint)");
    Count<int>::reset();

    cint2 = ccnt;
    Count<int>::dumpCounts("cint2=ccnt");
    Count<int>::reset();
    cout << "cint2=" << cint2 << endl;

    Array_<int,SmallIx> ismall0;
    Array_<int,SmallIx> ismall(3);
    Array_<int,SmallIx> imaxsz(data, data+4);
    cout << "ismall0=" << ismall0 << endl;
    cout << "ismall=" << ismall << endl;
    cout << "imaxsz=" << imaxsz << endl;

    new(ismall.raw_push_back()) int(27);
    cout << "ismall after raw_push_back():" << ismall << endl;

    cout << "sizeof(Array_<int,char>)=" << sizeof(Array_<int,char>) << endl;
    cout << "sizeof(Array_<int>)=" << sizeof(Array_<int>) << endl;
    cout << "sizeof(std::vector<int>)=" << sizeof(std::vector<int>) << endl;

    Array_<String, TestIx, 7> strings(6, "woohoo");
    cout << "strings=" << strings << endl;
    strings.push_back("last");
    cout << "strings=" << strings << endl;
    Array_<String, TestIx, 7>::const_reverse_iterator p = strings.rbegin();
    while (p != strings.rend())
        cout << " " << *p++;
    cout << endl;
}


int main() {
    SimTK_START_TEST("TestArray");

        SimTK_SUBTEST(testConstruction);

    SimTK_END_TEST();
}

