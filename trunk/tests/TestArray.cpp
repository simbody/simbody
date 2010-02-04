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
#include <sstream>
#include <iterator>
#include <iostream>
using std::cout;
using std::endl;
using std::cin;


using namespace SimTK;

// Output an std::vector<T>
template <class T>
std::ostream& operator<<(std::ostream& o, std::vector<T>& v) {
    o << '<';
    if (!v.empty()) {
        o << v.front();
        for (unsigned i=1; i < v.size(); ++i)
            o << ' ' << v[i];
    }
    return o << '>';
}

// Input an Array_<T>
template <class T, class X> inline std::istream&
operator>>(std::istream& i, Array_<T,X>& a) {
    if (a.isOwner()) {
        a.clear();
        while(!i.eof()) {
            a.push_back(); // default construct new element
            i >> a.back(); // input last element
        }
    } else { // non-owner
        typedef typename Array_<T,X>::size_type size_type;
        for (size_type e(0); e < a.size() && !i.eof(); ++e)
            i >> a[e];
    }
    return i;
} 
template <class T, class X> inline std::istream&
operator>>(std::istream& i, ArrayView_<T,X>& a) 
{   return i >> (Array_<T,X>&)a; }

template <class T>
class OtherArray_ : public Array_<T> {
public:
    typedef typename Array_<T>::size_type size_type;

    OtherArray_() : Array_<T>() {}
    OtherArray_(size_type n, const T& v) : Array_<T>(n,v) {}
};

SimTK_DEFINE_UNIQUE_INDEX_TYPE(TestIx);

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

namespace SimTK {
template <> struct NiceTypeName<SmallIx> {
    static const char* name() {return "SmallIx";}
};
}

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

// Instantiate the whole class to check for compilation problems.
template class Array_<int>;
template class Array_<std::string, unsigned char>;

#ifdef _MSC_VER // gcc 4.1.2 had trouble with this
// Instantiate templatized methods
typedef std::set<float>::const_iterator inputIt; // not a random access iterator

// Constructors.
template Array_<float,int>::Array_(const Array_<float,int>&);
template Array_<float,int>::Array_(const float*,const float*);

// Assignment.
template void 
Array_<float,int>::assign(const float*,const float*);
template void
Array_<double,int>::assign(const inputIt&, const inputIt&);
template Array_<float,int>& 
Array_<float,int>::operator=(const Array_<float,int>&);
template Array_<double,int>& 
Array_<double,int>::operator=(const std::vector<float>&);

// Insertion
template float*
Array_<float,int>::insert(float*, const float*, const float*);
template float*
Array_<float,short>::insert(float*, const inputIt&, const inputIt&);

// Comparison
template bool SimTK::operator==(const ArrayViewConst_<float,int>&, 
                                const ArrayViewConst_<float,unsigned>&);
#endif


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
    cout << "default constructed Array_<int> begin()=" << ismall0.begin()
         << " end()=" << ismall0.end() 
         << " capacity()=" << (int)ismall0.capacity() 
         << endl;

    std::vector<int> ivec0;
    cout << "default constructed std::vector<int>" 
         << " capacity()=" << ivec0.capacity() 
         << endl;

    Array_<int,SmallIx> ismall(3);
    Array_<int,SmallIx> imaxsz(data, data+4);
    cout << "ismall0=" << ismall0 << endl;
    cout << "ismall=" << ismall << endl;
    cout << "imaxsz=" << imaxsz << endl;

    new(ismall.raw_push_back()) int(27);
    cout << "ismall after raw_push_back():" << ismall << endl;

    SimTK_TEST_MUST_THROW_DEBUG(imaxsz.push_back()); // already full

    // Check null assignments.
    ismall = ismall0; // src is null
    ismall0 = imaxsz; // dest was null
    ismall = Array_<int,SmallIx>(); // both null

    cout << "sizeof(Array_<int,bool>)=" << sizeof(Array_<int,bool>) << endl;
    cout << "sizeof(Array_<int,char>)=" << sizeof(Array_<int,unsigned char>) << endl;
    cout << "sizeof(Array_<int,short>)=" << sizeof(Array_<int,unsigned short>) << endl;
    cout << "sizeof(Array_<int>)=" << sizeof(Array_<int>) << endl;
    cout << "sizeof(std::vector<int>)=" << sizeof(std::vector<int>) << endl;
    cout << "sizeof(Array_<int,long long>)=" << sizeof(Array_<int,long long>) << endl;

    Array_<String, TestIx> strings(6, "woohoo");
    cout << "strings=" << strings << endl;
    strings.push_back("last");
    for (int i=0; i<5; ++i) {
        strings.insert(strings.end(), 2, "ins" + String(i));
        cout << strings.size() << ":" << strings.capacity() 
                               << ":" << strings << endl;
    }
    cout << "strings=" << strings << endl;

    Array_<String, TestIx>::reverse_iterator p = strings.rbegin();
    while (p != strings.rend())
        cout << " " << *p++;
    cout << endl;

    const int ownerData[] = {7, 77, 777, 7777, 77777};
    std::vector<int> owner(ownerData, ownerData+5);
    std::vector<unsigned> unowner(owner.begin(), owner.end());
    Array_<int> shared; shared.shareData(&owner[1], &owner[4]);
    cout << "vector before=" << owner << endl;
    cout << "shared before=" << shared << endl;
    shared[2] = 29;
    cout << "shared after=" << shared << endl;
    cout << "vector after=" << owner << endl;
    cout << "shared(1,2)=" << shared(1,2) << endl;

    Array_<int> copyOfOwner(owner);
    cout << "copyOfOwner=" << copyOfOwner << endl;
    Array_<unsigned short,char> weirdCopy(owner);
    cout << "weirdCopy=" << weirdCopy << endl;
    copyOfOwner = unowner;
    cout << "copyOfOwner=unowner=" << copyOfOwner << endl;

    Array_<unsigned> shareOfUnowner(unowner, DontCopy());
    cout << "shareOfUnowner=" << shareOfUnowner << endl;

    shareOfUnowner(1,3) = Array_<unsigned>(3,88);
    cout << "shareOfUnowner=" << shareOfUnowner << endl;

    OtherArray_<int> oa(5, -4);
    cout << "oa=" << oa << endl;
}

static void toArray(const Array_<int>& a) {
    cout << "toArray=     "  << a  << "  &a[0]="  << &a[0] << endl;
}
static void toArrayView(const ArrayView_<int>& av) {
    cout << "toArrayView= "  << av  << " &av[0]="  << &av[0] << endl;
}
static void toArrayViewConst(const ArrayViewConst_<int>& ca) {
    cout << "toArrayViewConst="  << ca  << " &ca[0]="  << &ca[0] << endl;
}
void testConversion() {
    const int p[] = {1,2,3,4,5,6};
    std::vector<int> v(p,p+6);
    cout << "v=" << v << " &v[0]=" << &v[0] << endl;
    Array_<int> a(v);
    ArrayView_<int> av(v);
    ArrayViewConst_<int> ca(v);
    cout << "a= "  << a  << "  &a[0]="  << &a[0] << endl;
    cout << "av=" << av << " &av[0]=" << &av[0] << endl;
    cout << "ca=" << ca << " &ca[0]=" << &ca[0] << endl;

    toArray(ArrayView_<int>(v)); 
    toArrayView(v); 
    toArrayViewConst(v);
}

// ArrayView assignment can't change the size of the target ArrayView,
// and the semantics are elementwise assignment, not destruct-then-copy-
// construct as for resizeable Array (or std::vector) assignment.
void testArrayViewAssignment() {
    const int data[5] = {10, 100, -23, 4, -99};
    Array_<int> adata(data, data+5);    // copies of original data
    std::vector<int> vdata(data, data+5);

    Count<int>::reset();
    Array_< Count<int> > acnt(data, data+5);
    SimTK_TEST(!(  Count<int>::defCtor   ||Count<int>::copyCtor
                 ||Count<int>::copyAssign||Count<int>::dtor));
    SimTK_TEST(Count<int>::initCtor == 5);

    Count<int>::reset();
    acnt = adata; // clear() then construct from int
    SimTK_TEST(!(  Count<int>::defCtor   ||Count<int>::copyCtor
                 ||Count<int>::copyAssign));
    SimTK_TEST(Count<int>::dtor==5 && Count<int>::initCtor==5);

    Count<int>::reset();
    Array_< Count<int> > acopy(3); // default constructed
    SimTK_TEST(Count<int>::defCtor == 3);
    acopy = acnt; // destruct 3, copy construct 5
    SimTK_TEST(Count<int>::dtor==3 && Count<int>::copyCtor==5);

    Count<int>::reset();
    // this is an initialization, not an assignment
    ArrayView_< Count<int> > avcnt = acnt(1,2); // shares 2nd & 3rd elts
    SimTK_TEST(Count<int>::isReset()); // nothing should have happened
    SimTK_TEST(avcnt.size()==2);
    SimTK_TEST(avcnt[0]==acnt[1]&&avcnt[1]==acnt[2]);
    SimTK_TEST(&avcnt[0]==&acnt[1]&& &avcnt[1]==&acnt[2]);

    // This assignment should fail because the source has too many elements.
    SimTK_TEST_MUST_THROW(avcnt = adata);
    // This one should succeed, with 2 calls to Count<int>::op=(int).
    Count<int>::reset();
    avcnt.assign(adata.begin(), adata.begin()+2);
    SimTK_TEST(!(Count<int>::defCtor||Count<int>::copyCtor
                 ||Count<int>::copyAssign));
    SimTK_TEST(Count<int>::initAssign == 2);

    // This assignment should fail because of overlap between source and
    // destination.
    SimTK_TEST_MUST_THROW(avcnt = acnt(0,2));
    // But this succeeds because no overlap.
    Count<int>::reset();
    avcnt = acnt(3,2); // sets acnt(1,2)=acnt(3,2) with 2 copy assigns
    SimTK_TEST(!(Count<int>::copyCtor||Count<int>::dtor));
    SimTK_TEST(Count<int>::copyAssign == 2);
    //was: int data[5] = {10, 100, -23, 4, -99};
    int modified[5] = {10, 4, -99, 4, -99}; // should now be
    SimTK_TEST(avcnt[0]==4&&avcnt[1]==-99);
    SimTK_TEST(acnt == Array_<int>(modified, modified+5));

    // Check behavior of assign(first,last1) for input, forward, and
    // random access iterators.
    int someSpace[5] = {123, 1, 12, -9, 14};
    ArrayView_<int> avSpace(someSpace, someSpace+5);
    std::vector<int> aVec(avSpace.begin(), avSpace.end()); // copy
    std::set<int> aSet(avSpace.begin(), avSpace.end()); // copy & sort
    SimTK_TEST(&avSpace[0] == someSpace); //must be sharing space
    // Test fill first.
    avSpace = 19; SimTK_TEST(avSpace==Array_<int>(5, 19));
    avSpace.fill(-3); SimTK_TEST(avSpace==Array_<int>(5, -3));
    SimTK_TEST_MUST_THROW_DEBUG(avSpace.assign(12, 999));
    avSpace.assign(5,999); SimTK_TEST(avSpace==Array_<int>(5,999));
    // Assign from pointers
    avSpace.assign(data, data+5); SimTK_TEST(avSpace==vdata);
    avSpace = 999;
    // Assign from random_access_iterators
    avSpace.assign(vdata.begin(),vdata.end());SimTK_TEST(avSpace==vdata);
    avSpace = 999;
    // Assign from bidirectional_interator
    avSpace.assign(aSet.begin(), aSet.end());
    SimTK_TEST(avSpace==std::vector<int>(aSet.begin(),aSet.end()));

    SimTK_TEST_MUST_THROW(avSpace.assign(data, data+3));
    SimTK_TEST_MUST_THROW(avSpace.assign(vdata.begin(), vdata.begin()+3));
    std::set<int>::iterator sp = aSet.begin(); ++sp; ++sp;
    SimTK_TEST_MUST_THROW(avSpace.assign(aSet.begin(), sp));
}

// A bool index type is more or less useless in real life but was handy for 
// catching obscure implementation bugs having to do with index_type vs. 
// size_type.
void testBoolIndex() {
    SimTK_TEST((Array_<int,long>().empty()));
    SimTK_TEST(!(Array_<int,long>(1L,99).empty()));
    SimTK_TEST((Array_<int,long>(1L,99)[0] == 99));

    Array_<std::string, bool> wisdom(2);
    SimTK_TEST(wisdom[true] == ""); SimTK_TEST(wisdom[false] == "");
    wisdom[true]  = "this too shall pass";
    wisdom[false] = "don't worry it's not loaded";
    SimTK_TEST(wisdom.size() == 2);
    SimTK_TEST(wisdom.max_size() == 2);
    SimTK_TEST(wisdom.capacity() >= 2);
    SimTK_TEST(wisdom.allocated() == wisdom.capacity());
    SimTK_TEST(wisdom.data() != 0);
    SimTK_TEST(wisdom.data() == wisdom.cdata());
    SimTK_TEST(wisdom.begin() == wisdom.data());
    SimTK_TEST(wisdom.cbegin() == wisdom.cdata());
    SimTK_TEST(wisdom.end() == wisdom.begin()+2);
    SimTK_TEST(wisdom.cend() == wisdom.cbegin()+2);

    SimTK_TEST(wisdom[false] == "don't worry it's not loaded");
    SimTK_TEST(wisdom[true] == "this too shall pass");
    SimTK_TEST(wisdom.at(false) == "don't worry it's not loaded");
    SimTK_TEST(wisdom.at(true) == "this too shall pass");

    cout << "wisdom=" << wisdom << endl;
    cout << "wisdom(false,1)=" << wisdom(false,1) << endl;
    cout << "wisdom(true,1)=" << wisdom(true,1) << endl;
    cout << "wisdom(true,0)=" << wisdom(true,0) << endl;

    // Subarrays are fixed size; can't assign a 1-element vector to
    // a 2 element subarray.
    SimTK_TEST_MUST_THROW_DEBUG(
        wisdom(false,2) = std::vector<const char*>(1,"whatever"));

    const std::vector<const char*> vrel(2,"it's all relative");
    wisdom(false,2) = vrel;
    cout << "wisdom=" << wisdom << endl;

    // Test all the comparison operators Array vs. std::vector.
    SimTK_TEST(wisdom == vrel); SimTK_TEST(vrel == wisdom);
    SimTK_TEST(wisdom <= vrel); SimTK_TEST(vrel <= wisdom);
    SimTK_TEST(wisdom >= vrel); SimTK_TEST(vrel >= wisdom);
    SimTK_TEST(wisdom(false,1) < vrel);
    SimTK_TEST(wisdom(true,1) < vrel);
    SimTK_TEST(wisdom(0,0) < vrel);
    SimTK_TEST(wisdom(true,1) < vrel);
    SimTK_TEST(wisdom(false,1) != vrel);
    SimTK_TEST(wisdom != std::vector<const char*>(2,"it's all absolute"));
    wisdom[true] = "z comes after i";
    SimTK_TEST(wisdom > vrel); SimTK_TEST(vrel < wisdom);


    SimTK_TEST_MUST_THROW_DEBUG(wisdom(true,2));
    SimTK_TEST_MUST_THROW_DEBUG(wisdom.push_back("more brilliance"));
}

// It should be possible to assign to an Array_ from an std::set or std::map
// even though you can't subtract their bidirectional_iterators.
void testNonRandomIterator() {
    const int someInts[] = {30,40,10,20,30,7,5};
    std::set<int> iset(someInts, someInts+7);
    std::vector<int> sortUniq(iset.begin(), iset.end()); // the right answer

    Array_<int> iarr(iset.begin(), iset.end());
    iarr.assign(iset.begin(), iset.end()); // must increment to count
    SimTK_TEST(iarr == sortUniq);

    Array_<int> iarr2(sortUniq.begin(), sortUniq.end());
    iarr2.assign(sortUniq.begin(), sortUniq.end()); // can subtract iterators
    SimTK_TEST(iarr2 == sortUniq);
    iarr2.assign(iarr.begin(), iarr.end()); // can subtract pointers
    SimTK_TEST(iarr2 == sortUniq);

    // The standard requires this to match the constructor that creates
    // n copies of an initial value -- it must NOT match the templatized
    // InputIterator form because these are integral types.
    Array_<int> dummy1((char)3, 'A'); // 3*65
    SimTK_TEST(dummy1 == Array_<int>(3, (int)'A'));
    Array_<int> dummy2(4U, 129U);     // 4*129
    SimTK_TEST(dummy2 == Array_<int>(4, 129));

    // This should use the constant-time std::swap specialization that is 
    // provided in the Array.h header file.
    std::swap(dummy1, dummy2);
    SimTK_TEST(dummy2 == Array_<int>(3, (int)'A'));
    SimTK_TEST(dummy1 == Array_<int>(4, 129));

    // assign() and insert() should behave like the constructor.
    dummy1.assign((char)2, 'B');
    dummy1.insert(dummy1.begin()+1, (char)3, 'C');
    const int d1answer[] = {(int)'B',(int)'C',(int)'C',(int)'C',(int)'B'};
    SimTK_TEST((dummy1 == Array_<int,unsigned short>(d1answer, d1answer+5)));

    // Test fill().
    dummy1.fill(7);
    SimTK_TEST((dummy1 == Array_<int>(5, 7))); // i.e., 5 7's


    // This is too much data and should be detectable for any forward iterator.
    typedef Array_<int,SmallIx> AType;
    SimTK_TEST_MUST_THROW_DEBUG(
        AType small(iset.begin(), iset.end())); // bidirectional
    SimTK_TEST_MUST_THROW_DEBUG(
        AType small(sortUniq.begin(), sortUniq.end())); // random access
    SimTK_TEST_MUST_THROW_DEBUG(
        AType small(iarr.begin(), iarr.end())); // pointer

}

// Input iterators require special handling in the implementation because
// it can't be determined for them how many elements are in a range
// [first,last1) because to increment an iterator is to consume it. This is
// the only case where a bulk constructor, insert(), or assign() must be
// done with multiple space reallocations (basically like a series of 
// one-element push_back() calls).
void testInputIterator() {
    const int answerData[]={10,12,-14,5,203,-232,1,2,3,4};
    const Array_<int,char> answer(answerData,answerData+10);
    const Array_<int,short> smallAnswer(answerData+4, answerData+8);
    std::istringstream inp1("10 12 -14 5 203 -232 1 2 3 4");
    std::istringstream inp2("10 12 -14 5 203 -232 1 2 3 4");
    std::istringstream smallInp("203 -232 1 2"); // fits in SmallIx
    typedef std::istream_iterator<int> Iter;
    Iter p1(inp1), p2(inp2), psmall(smallInp);
    Array_<int> readin(p1, Iter()); // like begin(), end()
    SimTK_TEST(readin == answer);

    // This shouldn't work because there are too many elements.
    typedef Array_<int,SmallIx> SmallArray;
    SimTK_TEST_MUST_THROW_DEBUG(SmallArray tooSmall(p2, Iter()));

    // This should be OK.
    SmallArray okSmall(psmall, Iter());
    SimTK_TEST(okSmall == smallAnswer);

    Array_<float> farray;
    const float farray_ans1[] = {-1.5f,3e4f,.125f,11,4e-7f};
    std::istringstream fin1("-1.5 3e4 .125 11 4e-7");
    fin1 >> farray;
    SimTK_TEST(farray == std::vector<float>(farray_ans1, farray_ans1+5));

    // Replace middle three elements.
    ArrayView_<float> fmid(farray(1,3));
    const float farray_ans2[] = {-1.5f,910,920,9200,4e-7f};
    std::istringstream fin2("9.1e2 9.2e2 9.2e3");
    fin2 >> fmid;
    SimTK_TEST(farray == Array_<float>(farray_ans2, farray_ans2+5));

}

// Reduce the loop count by 50X in Debug.
static const int Outer = 500000
#ifndef NDEBUG
          / 50
#endif
                ;

static const int Inner = 1000;
void testSpeedStdVector() {
    std::vector<int> v; 
    v.reserve(Inner);

    for (int i=0; i < Outer; ++i) {
        v.clear();
        for (int i=0; i < Inner; ++i)
            v.push_back(i);
    }

    int sum;
    for (int i=0; i < Outer; ++i) {
        sum = i;
        for (unsigned i=0; i < v.size(); ++i)
            sum += v[i];
    }
    cout << "std::vector sum=" << sum << endl;
}

void testSpeedSimTKArray() {
    Array_<int> v; 
    v.reserve(Inner);

    for (int i=0; i < Outer; ++i) {
        v.clear();
        for (int i=0; i < Inner; ++i)
            v.push_back(i);
    }

    int sum;
    for (int i=0; i < Outer; ++i) {
        sum = i;
        for (int i=0; i < v.size(); ++i)
            sum += v[i];
    }
    cout << "Array sum=" << sum << endl;
}

void testNiceTypeName() {
    cout << "Is64BitPlatform=" << NiceTypeName<Is64BitPlatformType>::name() << endl;
    cout << "packed_size_type<bool>=" 
        << NiceTypeName<ArrayIndexPackType<bool>::packed_size_type>::name() << endl;
    cout << "packed_size_type<char>=" 
        << NiceTypeName<ArrayIndexPackType<char>::packed_size_type>::name() << endl;
    cout << "packed_size_type<signed char>=" 
        << NiceTypeName<ArrayIndexPackType<signed char>::packed_size_type>::name() << endl;
    cout << "packed_size_type<unsigned char>=" 
        << NiceTypeName<ArrayIndexPackType<unsigned char>::packed_size_type>::name() << endl;
    cout << "packed_size_type<short>=" 
        << NiceTypeName<ArrayIndexPackType<short>::packed_size_type>::name() << endl;
    cout << "packed_size_type<unsigned short>=" 
        << NiceTypeName<ArrayIndexPackType<unsigned short>::packed_size_type>::name() << endl;
    cout << "packed_size_type<int>=" 
        << NiceTypeName<ArrayIndexPackType<int>::packed_size_type>::name() << endl;
    cout << "packed_size_type<unsigned>=" 
        << NiceTypeName<ArrayIndexPackType<unsigned>::packed_size_type>::name() << endl;
    cout << "packed_size_type<long>=" 
        << NiceTypeName<ArrayIndexPackType<long>::packed_size_type>::name() << endl;
    cout << "packed_size_type<unsigned long long>=" 
        << NiceTypeName<ArrayIndexPackType<unsigned long long>::packed_size_type>::name() << endl;
    cout << NiceTypeName< Array_<String,char> >::name() << endl;
}

// The Array_ class is supposed to make better use of memory than does
// std::vector when the index type is smaller than a pointer.
void testMemoryFootprint() {
    // These conditions should apply on any 32- or 64-bit platform.
    SimTK_TEST(sizeof(Array_<int>)      <= sizeof(std::vector<int>));
    SimTK_TEST(sizeof(Array_<int,bool>) <  sizeof(std::vector<int>));
    SimTK_TEST(sizeof(Array_<int,char>) <  sizeof(std::vector<int>));
    SimTK_TEST(sizeof(Array_<int,signed char>)    <  sizeof(std::vector<int>));
    SimTK_TEST(sizeof(Array_<int,unsigned char>)  <  sizeof(std::vector<int>));
    SimTK_TEST(sizeof(Array_<int,short>)          <  sizeof(std::vector<int>));
    SimTK_TEST(sizeof(Array_<int,unsigned short>) <  sizeof(std::vector<int>));

    // Since an int is smaller than a pointer here we will do better than
    // any 3-pointer implementation. And we shouldn't be worse then normal
    // for long longs.
    if (Is64BitPlatform) {
        SimTK_TEST(sizeof(Array_<int,int>)       <  sizeof(std::vector<int>));
        SimTK_TEST(sizeof(Array_<int,unsigned>)  <  sizeof(std::vector<int>));
        SimTK_TEST(sizeof(Array_<int,long long>) <=  sizeof(std::vector<int>));
        SimTK_TEST(sizeof(Array_<int,unsigned long long>) <=  sizeof(std::vector<int>));
    }

    // We don't know if long will be 32 or 64 bit on any given 64 bit
    // implementation (it is 32 bits for MSVC and 64 for gcc). But it is
    // always 32 bits on a 32 bit implementation so we shouldn't be doing
    // any worse here.
    SimTK_TEST(sizeof(Array_<int,long>) <= sizeof(std::vector<int>));

    // Check that packing is working right.
    // ints and larger are treated the same for 32 vs 64. (longs are 
    // wobblers though so we don't check here)
    SimTK_TEST(sizeof(Array_<int>::packed_size_type)==sizeof(int));
    SimTK_TEST(sizeof(Array_<int,int>::packed_size_type)==sizeof(int));
    SimTK_TEST(sizeof(Array_<int,unsigned int>::packed_size_type)==sizeof(int));
    SimTK_TEST(sizeof(Array_<int,long long>::packed_size_type)==sizeof(long long));
    SimTK_TEST(sizeof(Array_<int,unsigned long long>::packed_size_type)==sizeof(long long));
    if (Is64BitPlatform) {
        // Small types are packed into an int on 64 bit platform.
        SimTK_TEST(sizeof(Array_<int,bool>::packed_size_type)==sizeof(int));
        SimTK_TEST(sizeof(Array_<int,char>::packed_size_type)==sizeof(int));
        SimTK_TEST(sizeof(Array_<int,signed char>::packed_size_type)==sizeof(int));
        SimTK_TEST(sizeof(Array_<int,unsigned char>::packed_size_type)==sizeof(int));
        SimTK_TEST(sizeof(Array_<int,short>::packed_size_type)==sizeof(int));
        SimTK_TEST(sizeof(Array_<int,unsigned short>::packed_size_type)==sizeof(int));
    } else { 
        // Small types are packed into a short on 32 bit platform.
        SimTK_TEST(sizeof(Array_<int,bool>::packed_size_type)==sizeof(short));
        SimTK_TEST(sizeof(Array_<int,char>::packed_size_type)==sizeof(short));
        SimTK_TEST(sizeof(Array_<int,signed char>::packed_size_type)==sizeof(short));
        SimTK_TEST(sizeof(Array_<int,unsigned char>::packed_size_type)==sizeof(short));
        SimTK_TEST(sizeof(Array_<int,short>::packed_size_type)==sizeof(short));
        SimTK_TEST(sizeof(Array_<int,unsigned short>::packed_size_type)==sizeof(short));
    }

    // Now we'll bravely insist that we know how these should be packed.
    if (Is64BitPlatform) {
        SimTK_TEST(sizeof(Array_<int>)==16);
        SimTK_TEST(sizeof(Array_<int,bool>)==16);
        SimTK_TEST(sizeof(Array_<int,char>)==16);
        SimTK_TEST(sizeof(Array_<int,signed char>)==16);
        SimTK_TEST(sizeof(Array_<int,unsigned char>)==16);
        SimTK_TEST(sizeof(Array_<int,short>)==16);
        SimTK_TEST(sizeof(Array_<int,unsigned short>)==16);
        SimTK_TEST(sizeof(Array_<int,int>)==16);
        SimTK_TEST(sizeof(Array_<int,unsigned>)==16);
        SimTK_TEST(sizeof(Array_<int,long>)<=24);
        SimTK_TEST(sizeof(Array_<int,unsigned long>)<=24);
        SimTK_TEST(sizeof(Array_<int,long long>)==24);
        SimTK_TEST(sizeof(Array_<int,unsigned long long>)==24);
    } else { // 32 bit platform
        SimTK_TEST(sizeof(Array_<int>)==12);
        SimTK_TEST(sizeof(Array_<int,bool>)==8);
        SimTK_TEST(sizeof(Array_<int,char>)==8);
        SimTK_TEST(sizeof(Array_<int,signed char>)==8);
        SimTK_TEST(sizeof(Array_<int,unsigned char>)==8);
        SimTK_TEST(sizeof(Array_<int,short>)==8);
        SimTK_TEST(sizeof(Array_<int,unsigned short>)==8);
        SimTK_TEST(sizeof(Array_<int,int>)==12);
        SimTK_TEST(sizeof(Array_<int,unsigned>)==12);
        SimTK_TEST(sizeof(Array_<int,long>)<=12);
        SimTK_TEST(sizeof(Array_<int,unsigned long>)<=12);
        // These don't make sense on a 32 bit platform, but they work. The
        // size will be 20 or 24 depending on how the compiler aligns the
        // 8-byte integers after the pointer.
        SimTK_TEST(sizeof(Array_<int,long long>)<=24);
        SimTK_TEST(sizeof(Array_<int,unsigned long long>)<=24);
    }
}

int main() {

    SimTK_START_TEST("TestArray");

        SimTK_SUBTEST(testArrayViewAssignment);
        SimTK_SUBTEST(testInputIterator);
        SimTK_SUBTEST(testNiceTypeName);
        SimTK_SUBTEST(testMemoryFootprint);
        SimTK_SUBTEST(testConstruction);
        SimTK_SUBTEST(testConversion);
        SimTK_SUBTEST(testBoolIndex);
        SimTK_SUBTEST(testNonRandomIterator);
        SimTK_SUBTEST(testSpeedStdVector);
        SimTK_SUBTEST(testSpeedSimTKArray);

    SimTK_END_TEST();
}

