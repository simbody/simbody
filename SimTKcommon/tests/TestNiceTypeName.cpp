/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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

/*
 * These are regression tests for the SimTK::NiceTypeName<T> class.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include <vector>
#include <sstream>
#include <iterator>
#include <iostream>
#include <utility>
#include <memory>
using std::cout;
using std::endl;
using std::cin;


using namespace SimTK;

template <class T>
class OtherArray_ : public Array_<T> {
public:
    typedef typename Array_<T>::size_type size_type;

    OtherArray_() : Array_<T>() {}
    OtherArray_(size_type n, const T& v) : Array_<T>(n,v) {}
};

SimTK_DEFINE_UNIQUE_INDEX_TYPE(TestIx);

// This index type has a max size of 4 for testing out-of-space
// checks.
class SmallIx {
public:
    SmallIx() : ix(0xff) {}
    explicit SmallIx(unsigned char i) : ix(i) {}

    SmallIx& operator++() 
    {   assert(ix<max_size()); ++ix; return *this;}
    SmallIx operator++(int) 
    {   assert(ix<max_size()); const SmallIx x=*this; ++ix; return x;}
    SmallIx& operator--() 
    {   assert(ix>0); --ix; return *this;}
    SmallIx operator--(int) 
    {   assert(ix>0); const SmallIx x=*this; ++ix; return x;}

    // These are required for any class to be used an index type.
    operator              unsigned char() const {return ix;}
    typedef unsigned char size_type;
    typedef signed char   difference_type;
    static size_type      max_size() {return 4;}
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

// This class is a T but augmented with counters that track the number
// of calls to constructors, assignment, and the destructor.
template <class T>
struct Count {
    Count() {++defCtor;}
    Count(const Count& c) : val(c.val) {++copyCtor;}
    Count& operator=(const Count& c) {val=c.val; ++copyAssign; return *this;}
    ~Count() {++dtor;}

    // Conversion from T.
    Count(const T& t) : val(t) {++initCtor;}
    // Assign from T.
    Count& operator=(const T& t) {val=t; ++initAssign; return *this;}

    enum Color {Red,Green,Blue};
    enum class Letter {A,B,C};


    bool operator==(const Count& other) const {return val==other.val;}
    bool operator!=(const Count& other) const {return val!=other.val;}

    static void dumpCounts(const char* msg) {
        cout << msg << ":";
        cout << " defCtor=" << Count<int>::defCtor;
        cout << " initCtor=" << Count<int>::initCtor;
        cout << " copyCtor=" << Count<int>::copyCtor;
        cout << " initAssign=" << Count<int>::initAssign;
        cout << " copyAssign=" << Count<int>::copyAssign;
        cout << " dtor=" << Count<int>::dtor;
        cout << endl;
    }

    static bool isReset() 
    {   return !(defCtor||initCtor||copyCtor||initAssign||copyAssign||dtor); }

    T val;

    static void reset() {defCtor=initCtor=copyCtor=initAssign=copyAssign=dtor=0;}
    static Counter defCtor;
    static Counter initCtor;
    static Counter copyCtor;
    static Counter initAssign;
    static Counter copyAssign;
    static Counter dtor;
};
template <class T> inline std::ostream&
operator<<(std::ostream& o, const Count<T>& c) {
    return o << c.val;
} 
template <class T> Counter Count<T>::defCtor;
template <class T> Counter Count<T>::initCtor;
template <class T> Counter Count<T>::copyCtor;
template <class T> Counter Count<T>::initAssign;
template <class T> Counter Count<T>::copyAssign;
template <class T> Counter Count<T>::dtor;

// Standalone tests of the method that is used by NiceTypeName<T>::namestr()
// to clean up the demangled names on various platforms.
void testCanonicalize() {
    // Standardize "unsigned int" to "unsigned"; get rid of extra spaces.
    SimTK_TEST(canonicalizeTypeName("class std :: vector < unsigned int >")
               == "std::vector<unsigned>");
    // OSX's stl like to throw in these extra namespaces.
    SimTK_TEST(canonicalizeTypeName("std:: __1 :: __23 :: set<T>")
               == "std::set<T>");

    // Should leaves spaces between words.
    SimTK_TEST(canonicalizeTypeName("lunch bucket")=="lunch bucket");
    // And keep funny looking namespaces if they aren't __digits.
    SimTK_TEST(canonicalizeTypeName("std::my__1::__23x::resigned char")
               == "std::my__1::__23x::resigned char");
}

void testBuiltins() {
    SimTK_TEST(NiceTypeName<bool>::namestr() == "bool");
    SimTK_TEST(NiceTypeName<signed char>::namestr() == "signed char");
    SimTK_TEST(NiceTypeName<unsigned char>::namestr() == "unsigned char");
    SimTK_TEST(NiceTypeName<short>::namestr() == "short");
    SimTK_TEST(NiceTypeName<unsigned short>::namestr() == "unsigned short");
    SimTK_TEST(NiceTypeName<int>::namestr() == "int");
    SimTK_TEST(NiceTypeName<unsigned int>::namestr() == "unsigned"); // abbr
    SimTK_TEST(NiceTypeName<unsigned>::namestr() == "unsigned");
    SimTK_TEST(NiceTypeName<long>::namestr() == "long");
    SimTK_TEST(NiceTypeName<unsigned long>::namestr() == "unsigned long");
    SimTK_TEST(NiceTypeName<long long>::namestr() == "long long");
    SimTK_TEST(NiceTypeName<unsigned long long>::namestr() 
               == "unsigned long long");
    SimTK_TEST(NiceTypeName<float>::namestr() == "float");
    SimTK_TEST(NiceTypeName<double>::namestr() == "double");
    SimTK_TEST(NiceTypeName<long double>::namestr() == "long double");
    SimTK_TEST(NiceTypeName<std::complex<float>>::namestr() 
               == "std::complex<float>");
    SimTK_TEST(NiceTypeName<std::complex<double>>::namestr() 
               == "std::complex<double>");
    SimTK_TEST(NiceTypeName<std::complex<long double>>::namestr() 
               == "std::complex<long double>");

    // xmlstr should be the same as namestr for non-templatized types
    SimTK_TEST(NiceTypeName<bool>::xmlstr() == NiceTypeName<bool>::namestr()); 
    SimTK_TEST(NiceTypeName<unsigned int>::xmlstr() 
               == NiceTypeName<unsigned int>::namestr());
    SimTK_TEST(NiceTypeName<long long>::xmlstr() 
               == NiceTypeName<long long>::namestr());

    // xmlstr should replace the brackets for templatized types
    SimTK_TEST(NiceTypeName<std::complex<float>>::xmlstr() 
               == "std::complex{float}");
    SimTK_TEST(NiceTypeName<std::complex<double>>::xmlstr() 
               == "std::complex{double}");
    SimTK_TEST(NiceTypeName<std::complex<long double>>::xmlstr() 
               == "std::complex{long double}");

}

namespace LocalNS {
enum MyEnum {One,Two};
enum class YourEnumClass {Three,Four};
}

void testEnums() {
    using namespace LocalNS;

    SimTK_TEST(NiceTypeName<MyEnum>::namestr() == "LocalNS::MyEnum");
    SimTK_TEST(NiceTypeName<YourEnumClass>::namestr() 
               == "LocalNS::YourEnumClass");

    SimTK_TEST(NiceTypeName<Count<double>>::namestr() == "Count<double>");
    SimTK_TEST(NiceTypeName<Count<double>::Color>::namestr() 
               == "Count<double>::Color");
    SimTK_TEST(NiceTypeName<Count<double>::Letter>::namestr() 
               == "Count<double>::Letter");
}

// Custom NiceTypeName for SmallIx: CustomSmallIxName
namespace SimTK {
template <> struct NiceTypeName<SmallIx> {
    static const char* name() {return "CustomSmallIxName";}
    static const std::string& namestr() 
    {   static const std::string ns(name()); return ns; }
    static const std::string& xmlstr() {return namestr();}
};

template <class T> 
struct NiceTypeName<OtherArray_<T>> {
    static const char* name() {return typeid(OtherArray_<T>).name();}
    static const std::string& namestr() 
    {   static const std::string ns
            ("OtherArray_<" + NiceTypeName<T>::namestr() + ">");
        return ns; }
    static const std::string& xmlstr() 
    {   static const std::string xs = encodeTypeNameForXML(namestr());
        return xs; }
};
}

void testArrayNames() {
    SimTK_TEST((NiceTypeName<Array_<String,char>>::namestr()
                == "SimTK::Array_<SimTK::String,char>"));
    SimTK_TEST((NiceTypeName<Array_<String,char>>::xmlstr()
                == "SimTK::Array_{SimTK::String,char}"));
    SimTK_TEST((NiceTypeName<ArrayView_<int>>::namestr()
                == "SimTK::ArrayView_<int,unsigned>"));
    SimTK_TEST((NiceTypeName<ArrayViewConst_<char>>::namestr()
                == "SimTK::ArrayViewConst_<char,unsigned>"));
    SimTK_TEST((NiceTypeName<Array_< Count<int> > >::namestr()
                == "SimTK::Array_<Count<int>,unsigned>"));
    SimTK_TEST((NiceTypeName<Array_< Count<int> > >::xmlstr()
                == "SimTK::Array_{Count{int},unsigned}"));
    SimTK_TEST((NiceTypeName<SmallIx>::namestr()
                == "CustomSmallIxName")); 
    SimTK_TEST((NiceTypeName<OtherArray_<int>>::namestr()
                == "OtherArray_<int>"));
    SimTK_TEST((NiceTypeName<OtherArray_<SmallIx>>::namestr()
                == "OtherArray_<CustomSmallIxName>"));
    SimTK_TEST(NiceTypeName<TestIx>::namestr() == "TestIx");
}

void testSTLNames() {
    SimTK_TEST(NiceTypeName<std::string>::namestr()
                == "std::string");
    SimTK_TEST((NiceTypeName<std::vector<int>>::namestr()
                == "std::vector<int,std::allocator<int>>"));
}


namespace {
// This is implicitly convertible to TextIx.
class SubTestIx : public TestIx {
public: explicit SubTestIx(int ix) : TestIx(ix) {}
};
}

// Output for anonymous or function-local type won't necessarily be platform 
// independent (and may be ugly), so we'll just write them out here for 
// inspection.
void testAnonymous() {
    cout << "These won't be nice:\n";
    cout << "(1) anonymous::SubTestIx gives " 
         << NiceTypeName<SubTestIx>::namestr() << endl;

    struct MyLocalType {
        int i;
    };

    cout << "(2) function local::MyLocalType gives " 
         << NiceTypeName<MyLocalType>::namestr() << endl;
}

int main() {

    SimTK_START_TEST("TestArray");

        SimTK_SUBTEST(testCanonicalize);
        SimTK_SUBTEST(testBuiltins);
        SimTK_SUBTEST(testEnums);
        SimTK_SUBTEST(testArrayNames);
        SimTK_SUBTEST(testSTLNames);
        SimTK_SUBTEST(testAnonymous);

    SimTK_END_TEST();
}

