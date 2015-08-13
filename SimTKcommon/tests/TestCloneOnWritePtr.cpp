/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

/* Test for proper functioning of the copy-on-write smart pointer class
SimTK::CloneOnWritePtr, and the null-on-copy variant of std::unique_ptr 
called SimTK::NullOnCopyUniquePtr. */

#include "SimTKcommon.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <utility>
#include <memory>
#include <vector>

using namespace SimTK;
using std::cout; using std::endl; using std::string; using std::unique_ptr;

// C++14
// using namespace std::literals;

class Base {
public:
    explicit Base(const std::string& n) : m_name(n) {
        ++m_constructions;
        //cout << "Base(" << getName() << ") @" << this << endl;
    }
    Base(const Base& src) : m_name(src.m_name) {
        ++m_copies;
        //cout << "CopyCtor Base(" << src.getName() << ") @" << &src 
        //     << " --> @" << this << endl;
    }
    virtual ~Base() {
        ++m_destructions;
        //cout << "~Base(" << getName() << ") @" << this << endl;
    }

    virtual Base* clone() const = 0;
    virtual int getValue() const = 0;
    virtual int& updValue() = 0;
    void setValue(int v) {updValue() = v;}
    const char* getName() const {return m_name.c_str();}

    static void dumpStats(const std::string& msg) {
        printf("dumpStats(%s):\n", msg.c_str());
        printf("# constructions=%d\n", m_constructions);
        printf("# copies=%d\n", m_copies);
        printf("# destructions=%d\n", m_destructions);
        printf("# alive: %d\n", getNumAlive());
    }

    // Return the net number of Base objects that have been constructed
    // but not yet destructed.
    static int getNumAlive() {
        return m_constructions+m_copies - m_destructions;
    }

    static int m_constructions;
    static int m_destructions;
    static int m_copies;
private:
    std::string m_name;

};
int Base::m_constructions = 0;
int Base::m_destructions = 0;
int Base::m_copies = 0;

std::ostream& operator<<(std::ostream& o, const Base& b) {
    o << "Base: " << b.getName() << "=" << b.getValue();
    return o;
}

template <class T>
std::ostream& operator<<(std::ostream& o,
                         const CloneOnWritePtr<T>& p) {
    o << "CloneOnWritePtr p=" << p.get() 
      << " use_count=" << p.use_count() << endl;
    if (p.empty()) o << "  EMPTY" << endl;
    else o << "  Obj: " << *p << endl; 
    return o;
}

class Derived1 : public Base {
public:
    explicit Derived1(const std::string& n) : Base(n) {
        //printf("Derived1() %s\n", getName());
    }
    ~Derived1() {/*printf("~Derived1(%s)\n", getName());*/}
private:
};

class Derived2 : public Base {
public:
    explicit Derived2(const std::string& n, int v) : Base(n), m_val2(v) {
        //printf("Derived2() %s\n", getName());
    }
    ~Derived2() {/*printf("~Derived2(%s)\n", getName());*/}

    Derived2* clone() const override {
        //printf("Derived2::clone(%llx) %s\n", 
        //       (unsigned long long)this, getName());
        return new Derived2(*this);
    }
    int getValue() const override {return m_val2;}
    int& updValue() override {return m_val2;}
private:
    int m_val2;
};

class Sub1 : public Derived1 {
public:
    explicit Sub1(const std::string& n, int v) : Derived1(n), m_val1(v)  {
        //printf("Sub1(%s,%d)\n", getName(), m_val1);
    }
    ~Sub1() {/*printf("~Sub1(%s)\n", getName());*/}
    Sub1* clone() const override {
        //printf("Sub1::clone(%llx) %s\n", 
        //       (unsigned long long)this, getName());
        return new Sub1(*this);
    }
    int getValue() const override {return m_val1;}
    int& updValue() override {return m_val1;}
private:
    int m_val1;
};

void testEmpty() {
    CloneOnWritePtr<Base> p, pp;
    SimTK_TEST(p.empty() && !p.unique());
    SimTK_TEST(p.use_count()==0);
    SimTK_TEST(!p);
    SimTK_TEST(!p.get() && !p.upd());
    SimTK_TEST(!p.release());
    SimTK_TEST_MUST_THROW_DEBUG(p.getRef());
    SimTK_TEST_MUST_THROW_DEBUG(p.updRef());
    SimTK_TEST_MUST_THROW_DEBUG(p->getValue());
    SimTK_TEST_MUST_THROW_DEBUG((*p).getValue());
    p.detach(); // shouldn't do anything
    p.swap(pp);
    std::swap(p,pp);
    SimTK_TEST(p.empty() && pp.empty());

    CloneOnWritePtr<Sub1> q(nullptr);
    SimTK_TEST(q.empty() && !q.unique() && q.use_count()==0);

    Derived2* dp = 0;
    CloneOnWritePtr<Derived2> d2(dp);
    SimTK_TEST(d2.empty() && !d2.unique() && d2.use_count()==0);

    // Check relational operators applied to empty Ptr.
    SimTK_TEST(!d2);
    SimTK_TEST(d2==nullptr && nullptr==d2);
    SimTK_TEST(d2==p && p==q); // these are all empty
    SimTK_TEST(!(d2!=nullptr || nullptr!=d2));
    SimTK_TEST(d2 >= nullptr && d2 <= nullptr && d2 >= p && d2 <= p);
    SimTK_TEST(!(d2 > nullptr || d2 < nullptr || d2 < p || d2 > p));

    p = pp; // copy assign
    p = d2; // copy assign compatible type
    p = std::move(pp); // move assign
    p = std::move(d2); // move assign compatible type
    SimTK_TEST(p.empty());
}

void testAllocate() {
    SimTK_TEST(Base::getNumAlive() == 0);

    CloneOnWritePtr<Base> p(new Sub1("first", 1));
    CloneOnWritePtr<Base> q(new Sub1("second", 2));
    CloneOnWritePtr<Base> r(new Derived2("d2", 999));
    CloneOnWritePtr<Base> e;

    SimTK_TEST(Base::getNumAlive() == 3);

    SimTK_TEST(p.unique() && !p.empty() && p.use_count()==1);
    SimTK_TEST(q.unique() && !q.empty() && q.use_count()==1);
    SimTK_TEST(r.unique() && !r.empty() && r.use_count()==1);
    SimTK_TEST(p != q && p != r && q != r);
    SimTK_TEST(e.empty());

    // Some stack-allocated objects.
    Sub1 localSub1("localSub1", 111);
    Derived2 localDerived2("localDerived2", 222);

    SimTK_TEST(Base::getNumAlive() == 5);

    CloneOnWritePtr<Sub1> psub1(localSub1);
    e = localDerived2;

    SimTK_TEST(psub1.unique());
    SimTK_TEST(e.unique());
    SimTK_TEST(Base::getNumAlive() == 7);


    // This should destroy p's holding, then share with psub1.
    p = psub1;
    SimTK_TEST(Base::getNumAlive() == 6);
    SimTK_TEST(!p.unique() && !psub1.unique());
    SimTK_TEST(p.use_count()==2 && psub1.use_count()==2);
    SimTK_TEST(p == psub1 && p != e);

    // Take over ownership of psub1's pointer. Since that is currently
    // shared with p, this should first create a new copy.
    Base* raw = psub1.release();
    SimTK_TEST(Base::getNumAlive() == 7);
    SimTK_TEST(psub1.empty() && p.unique() && raw);

    // Create a new Ptr to take over responsibility for takepsub1.
    // No allocation should occur.
    CloneOnWritePtr<Base> takeover(raw); raw=nullptr;
    SimTK_TEST(Base::getNumAlive() == 7);
    SimTK_TEST(takeover.unique());

    // Take the pointer back and then put it back in using assignment.
    raw = takeover.release();
    SimTK_TEST(Base::getNumAlive() == 7); // no allocation
    takeover = raw; raw=nullptr;
    SimTK_TEST(Base::getNumAlive() == 7); // still no allocation
    SimTK_TEST(takeover.unique());

    // This should invoke the move constructor, so only one allocation should
    // occur.
    CloneOnWritePtr<Base> cowd2 = CloneOnWritePtr<Derived2>(localDerived2);
    SimTK_TEST(Base::getNumAlive() == 8); // still no allocation
    cowd2 = nullptr; // should be same as reset()
    SimTK_TEST(Base::getNumAlive()==7 && cowd2.empty());

    CloneOnWritePtr<Derived2> d2Ptr(localDerived2);
    SimTK_TEST(Base::getNumAlive()==8 && d2Ptr.unique());

    // This should invoke the copy constructor, which is shared so this
    // should not require any allocations.
    CloneOnWritePtr<Base> cowd2again = d2Ptr;
    SimTK_TEST(Base::getNumAlive()==8);
    SimTK_TEST(cowd2again.use_count()==2 && d2Ptr.use_count()==2);
    SimTK_TEST(cowd2again == d2Ptr);

    // Writing to this should cause it to detach.
    cowd2again->setValue(3);
    SimTK_TEST(Base::getNumAlive()==9 && cowd2again.unique() && d2Ptr.unique());
    // The -> op could cause detach here but shouldn't because use_count==1.
    SimTK_TEST(cowd2again->getValue()==3 && Base::getNumAlive()==9);

    p=q=r; 
    SimTK_TEST(Base::getNumAlive()==7 && p.use_count()==3);
    SimTK_TEST(p==q && p==r && q==r);

    // This shouldn't detach since we're getting the const pointer.
    int rval = r.get()->getValue();
    SimTK_TEST(Base::getNumAlive()==7 && rval==999);

    // This *will* detach (unfortunately) since r is non-const.
    rval = r->getValue();
    SimTK_TEST(Base::getNumAlive()==8 && rval==999 && r.unique());
    SimTK_TEST(p.use_count()==2 && q.use_count()==2 && p==q && p!=r);

    // This should invoke const -> so should not detach.
    int qval = static_cast<const decltype(q)>(q)->getValue();
    SimTK_TEST(Base::getNumAlive()==8 && qval==999 && q.use_count()==2);

    (*q).updValue() = 111; // Should detach.
    qval = q.getRef().getValue();
    SimTK_TEST(Base::getNumAlive()==9 && qval==111 && q.unique());
    SimTK_TEST(r.get()->getValue()==999 && p.get()->getValue()==999 
               && r.unique() && p.unique() && q.unique());

    CloneOnWritePtr<Base> bptr; //empty

    bptr = r; // copy assign; no allocation
    SimTK_TEST(Base::getNumAlive()==9 && bptr.use_count()==2 && bptr==r);

    bptr.reset(); // bp is empty again; nothing deallocated
    SimTK_TEST(Base::getNumAlive()==9 && bptr.empty() && r.unique());

    bptr = std::move(r); // move assignment
    SimTK_TEST(Base::getNumAlive()==9 && bptr.unique() && r.empty());

    CloneOnWritePtr<Base> bptr2 = bptr; // copy construction
    SimTK_TEST(Base::getNumAlive()==9 && bptr2.use_count()==2 && bptr2==bptr);
    bptr2.detach(); // separate bptr2 and bptr
    SimTK_TEST(Base::getNumAlive()==10 && bptr2.unique() && bptr.unique());

    bptr2->setValue(101); bptr->setValue(-102);
    SimTK_TEST(Base::getNumAlive()==10
               && bptr2->getValue()==101 && bptr->getValue()==-102);
    std::swap(bptr, bptr2);
    SimTK_TEST(Base::getNumAlive()==10
               && bptr2->getValue()==-102 && bptr->getValue()==101);

    bptr=bptr2; // make them share again
    SimTK_TEST(Base::getNumAlive()==9);

    bptr2.reset(new Sub1("devilish", 666));
    SimTK_TEST(Base::getNumAlive()==10 && bptr2.get()->getValue()==666);
    SimTK_TEST(bptr.get()->getValue()==-102 && bptr.unique() && bptr2.unique());

    // Make sure upd() causes a detach.
    bptr=bptr2; // back to sharing
    SimTK_TEST(Base::getNumAlive()==9);

    bptr.upd()->setValue(999);
    SimTK_TEST(Base::getNumAlive()==10 && bptr.unique() && bptr2.unique());
    SimTK_TEST(bptr.get()->getValue()==999 && bptr2.get()->getValue()==666);

    //Test self-assignment and self move.
    bptr=bptr2; // sharing
    SimTK_TEST(Base::getNumAlive()==9);

    bptr=bptr; // should do nothing
    SimTK_TEST(Base::getNumAlive()==9 && bptr.use_count()==2);

    bptr=bptr2; // already sharing; should do nothing
    SimTK_TEST(Base::getNumAlive()==9 && bptr.use_count()==2);

    bptr=std::move(bptr2); // should reduce use count, empty bptr2
    SimTK_TEST(Base::getNumAlive()==9 && bptr.unique() && bptr2.empty());

    bptr=std::move(bptr); // self move should do nothing
    SimTK_TEST(Base::getNumAlive()==9 && bptr.unique());

}

// Call this at the end after all the destructors should have been
// called for anything allocated in the other tests. The Base class
// has been counting them.
void testForLeaks() {
    Base::dumpStats("FINAL");
    SimTK_TEST(Base::getNumAlive() == 0);
}

//#define TRACE(fmt, ...) printf(fmt, __VA_ARGS__)
#define TRACE(fmt, ...)
class Thing {
public:
    explicit Thing(int i) : t(i) {TRACE("construct T @0x%llx\n", 
                                      (unsigned long long)this);}
    ~Thing() { TRACE("destruct T @0x%llx\n", (unsigned long long)this); }
    Thing(const Thing& s) : t(s.t) { TRACE("copy construct T\n"); }
    Thing(Thing&& s) : t(std::move(s.t)) { TRACE("move construct T\n"); }
    Thing& operator=(const Thing& s) {
        t = s.t; TRACE("copy assign T\n");
        return *this;
    }
    Thing& operator=(Thing&& s) {
        t = s.t; TRACE("move assign T\n");
        return *this;
    }
    int getVal() const {return t;}
private:
    int t;
};

std::ostream& operator<<(std::ostream& o, const Thing& t) {
    return o << t.getVal();
}

void testNullOnCopyUniquePtr() {
    std::unique_ptr<Thing> u(new Thing(-3));
    NullOnCopyUniquePtr<Thing> pu(std::move(u));
    SimTK_TEST(pu->getVal() == -3 && !u);
    u = std::move(pu); // should move it back
    SimTK_TEST(u->getVal() == -3 && !pu);
    NullOnCopyUniquePtr<Thing> au;
    au = std::move(u); 
    SimTK_TEST(au->getVal() == -3 && !u);

    // Test swap.
    NullOnCopyUniquePtr<Thing> p(new Thing(3));
    NullOnCopyUniquePtr<Thing> q(new Thing(33));
    SimTK_TEST(p->getVal()==3 && q->getVal()==33);
    const Thing* pp = p.get(); const Thing* qp = q.get();
    swap(p,q);
    SimTK_TEST(q->getVal()==3 && p->getVal()==33);
    SimTK_TEST(q.get() == pp && p.get() == qp);

    NullOnCopyUniquePtr<Thing> np;
    NullOnCopyUniquePtr<Thing> npn(nullptr);
    SimTK_TEST(!np && !npn && np==nullptr && npn==nullptr);
    SimTK_TEST(nullptr==np && nullptr==npn);
    SimTK_TEST(np <= nullptr && np >= nullptr);
    SimTK_TEST(!(np > nullptr || np < nullptr));
    SimTK_TEST(np != q);

    NullOnCopyUniquePtr<Thing> cp = p; // copy constructor
    SimTK_TEST(!cp && p.get() == qp);

    NullOnCopyUniquePtr<Thing> ap;
    ap = p; // copy assignment
    SimTK_TEST(!ap && p.get() == qp);

    NullOnCopyUniquePtr<Thing> mp = std::move(p); // move constructor
    SimTK_TEST(mp->getVal()==33 && mp.get()==qp && !p);
    ap = std::move(mp); // move assignment
    SimTK_TEST(ap->getVal()==33 && ap.get()==qp && !mp);

    ap = nullptr; // null assignment inherited from std::unique_ptr
    SimTK_TEST(!ap);

    // Hash should be the same for pointer, unique_ptr, and our class.
    size_t hval = std::hash<Thing*>()(q.get());
    size_t hvalu = std::hash<std::unique_ptr<Thing>>()
                                    (static_cast<std::unique_ptr<Thing>&>(q));
    SimTK_TEST(hvalu == hval);

    size_t hvaln = std::hash<NullOnCopyUniquePtr<Thing>>()(q);
    SimTK_TEST(hvaln == hvalu);

    std::hash<NullOnCopyUniquePtr<Thing>> myhash;
    size_t hvaln2 = myhash(q);
    SimTK_TEST(hvaln2 == hvaln);
}

class UsesResetOnCopy {
public:
    ResetOnCopy<int>                    defint;
    ResetOnCopy<char>                   charZ = 'z';
    ResetOnCopy<string>                 defstr;
    ResetOnCopy<string>                 strHelloC = "hello";    // char* literal
    //ResetOnCopy<string>                 strGoodbyeS = "goodbye"s; // string literal C++14
    ResetOnCopy<short>                  shArr[3] = {9,8,7};
    ResetOnCopy<SubsystemIndex>         subIx{5};
    ResetOnCopy<std::vector<string>>    vstr {"one", "two", "three"};
    ResetOnCopy<unique_ptr<Array_<int>>> up{new Array_<int>({1,2,3})};

    void checkHasInitialValues() const {
        SimTK_TEST(defint == 0 && charZ == 'z');
        SimTK_TEST(defstr.empty());
        SimTK_TEST(strHelloC == "hello");
        //SimTK_TEST(strGoodbyeS == "goodbye");
        SimTK_TEST(shArr[0]==9 && shArr[1]==8 && shArr[2]==7);
        SimTK_TEST(subIx == 5);
        SimTK_TEST(vstr.size() == 3);
        SimTK_TEST(vstr[0]=="one" && vstr[2]=="three");
        SimTK_TEST((*up).size() == 3);
        SimTK_TEST((*up)[0]==1 && (*up)[2]==3);
    }

    void checkHasResetValues() const {
        SimTK_TEST(defint == 0 && charZ == '\0');
        SimTK_TEST(defstr.empty());
        SimTK_TEST(strHelloC.empty());
        //SimTK_TEST(strGoodbyeS.empty());
        SimTK_TEST(shArr[0]==0 && shArr[1]==0 && shArr[2]==0);
        SimTK_TEST(!subIx.isValid());
        SimTK_TEST(vstr.empty());
        SimTK_TEST(!up);
    }
};

void testResetOnCopy() {
    UsesResetOnCopy r1;
    r1.checkHasInitialValues();

    UsesResetOnCopy r2(r1); // copy construction
    r2.checkHasResetValues();

    r2 = r1; // copy assignment
    r2.checkHasResetValues();

    r2 = std::move(r1); // move assignment (r1 is damaged now)
    r2.checkHasInitialValues();

    UsesResetOnCopy r3(std::move(r2)); // move construction
    r3.checkHasInitialValues();

    string str1("hi there");
    ResetOnCopy<string> rcstr2(str1);
    SimTK_TEST(rcstr2 == "hi there");

    ResetOnCopy<int> rint(3);
    SimTK_TEST(rint == 3);

    ResetOnCopy<const double*> dstar;
    SimTK_TEST(dstar == nullptr);

    ResetOnCopy<double> rdub;
    SimTK_TEST(rdub == 0.);

    ResetOnCopy<string> mystr("hello");
    SimTK_TEST(mystr == "hello");

    mystr = std::move(mystr); // self-move assignment should do nothing
    SimTK_TEST(mystr == "hello");

    mystr = mystr; // copy assignment should clear
    SimTK_TEST(mystr.empty());

    ResetOnCopy<unique_ptr<double>> updub;
    SimTK_TEST(updub.get() == nullptr);

    updub.reset(new double(1.25));
    SimTK_TEST(updub && *updub == 1.25);

    // This is a move assignment since the RHS is a temporary rvalue. Copy
    // assignment wouldn't compile;
    double* dp = new double(3.125); 
    updub = unique_ptr<double>(dp); // should transfer the heap object
    SimTK_TEST(updub.get()==dp && *updub == 3.125);

    // This is copy construction so should default initialize and leave source
    // untouched.
    ResetOnCopy<unique_ptr<double>> updub2(updub);
    SimTK_TEST(!updub2 && updub.get()==dp && *updub == 3.125);

    // Move construction should move the heap object and leave the source
    // empty.
    ResetOnCopy<unique_ptr<double>> updub3(std::move(updub));
    SimTK_TEST(updub3.get()==dp && *updub3 == 3.125);
    SimTK_TEST(!updub);

    ResetOnCopy<SubsystemIndex> mysix(6);
    SimTK_TEST(mysix == 6);

    ResetOnCopy<SubsystemIndex> newsix(mysix); // copy construct
    SimTK_TEST(!newsix.isValid());

    newsix = SubsystemIndex(7); // ordinary assignment to contained type
    SimTK_TEST(newsix == 7);

    newsix = mysix; // copy assignment; should reinitialize
    SimTK_TEST(!newsix.isValid());

    newsix = std::move(mysix); // move assignment
    SimTK_TEST(newsix == 6);

    ResetOnCopy<int> arr[3];
    SimTK_TEST(arr[0]==0 && arr[1]==0 && arr[2]==0);
    ResetOnCopy<double> arr2[3] {9.25,8,7};
    SimTK_TEST(arr2[0]==9.25 && arr2[1]==8 && arr2[2]==7);
}

class UsesReinitOnCopy {
public:
    enum Color {Red, Green, Blue};

    ReinitOnCopy<char>                   charZ = 'z';
    ReinitOnCopy<string>                 strHelloC = "hello";    // char* literal
    // ReinitOnCopy<string>                 strGoodbyeS = "goodbye"s; // string literal C++14
    ReinitOnCopy<short>                  shArr[3] = {9,8,7};
    ReinitOnCopy<SubsystemIndex>         subIx{5};
    ReinitOnCopy<std::vector<string>>    vstr {"one", "two", "three"};
    ReinitOnCopy<Color>                  color{Blue};

    void checkHasInitialValues() const {
        SimTK_TEST(charZ == 'z');
        SimTK_TEST(strHelloC == "hello");
        //SimTK_TEST(strGoodbyeS == "goodbye");
        SimTK_TEST(shArr[0]==9 && shArr[1]==8 && shArr[2]==7);
        SimTK_TEST(subIx == 5);
        SimTK_TEST(vstr.size() == 3);
        SimTK_TEST(vstr[0]=="one" && vstr[2]=="three");
        SimTK_TEST(color == Blue);
    }

    void changeAll() {
        charZ = 'y';
        strHelloC = "something";
        //strGoodbyeS = "a string"s;
        // shArr[0] = shArr[1] = shArr[2] = -123 won't work: two of those are
        // copy assignments so get reinitialized instead.
        shArr[0] = -123; shArr[1] = -123; shArr[2] = -123;
        subIx.invalidate();
        vstr = std::vector<string>{"breakfast", "lunch"};
        color = Red;
    }

    void checkHasChanged() const {
        SimTK_TEST(charZ == 'y');
        SimTK_TEST(strHelloC == "something");
        //SimTK_TEST(strGoodbyeS == "a string");
        SimTK_TEST(shArr[0]==-123 && shArr[1]==-123 && shArr[2]==-123);
        SimTK_TEST(!subIx.isValid());
        SimTK_TEST(vstr.size() == 2);
        SimTK_TEST(vstr[0]=="breakfast" && vstr[1]=="lunch");
        SimTK_TEST(color == Red);
    }

};

void testReinitOnCopy() {
    UsesReinitOnCopy r1;
    r1.checkHasInitialValues();
    
    // Change every field and make sure it changes back on copy.
    r1.changeAll();
    r1.checkHasChanged();

    UsesReinitOnCopy r2(r1); // copy construction; should reinit
    r2.checkHasInitialValues();

    r2.changeAll();
    r2.checkHasChanged();

    r2 = r1; // copy assignment; should reinit.
    r2.checkHasInitialValues();

    r2 = std::move(r1); // move assignment; r2 gets changed values
    r2.checkHasChanged();

    UsesReinitOnCopy r3(std::move(r2)); // move construction; r3 gets changed
    r3.checkHasChanged();


    enum Color {Red=1,Green,Blue}; // 0 isn't one of the enumerators
    ReinitOnCopy<Color> myColor{Blue};
    SimTK_TEST(myColor == Blue);
    myColor = Green;
    SimTK_TEST(myColor == Green);

    ReinitOnCopy<Color> myColor2(myColor); // copy construction; reinit
    SimTK_TEST(myColor2 == Blue);

    ReinitOnCopy<double> d{123.};
    SimTK_TEST(d == 123.);
    d = 3.125;
    SimTK_TEST(d == 3.125);

    ReinitOnCopy<double> dd(d); // copy construction; reinit
    SimTK_TEST(dd == 123.);

    ReinitOnCopy<string> mystr1 = "unknown";
    SimTK_TEST(mystr1 == "unknown" && mystr1.getReinitValue() == "unknown");
    mystr1 = "now we know";
    SimTK_TEST(mystr1 == "now we know" && mystr1.getReinitValue() == "unknown");

    ReinitOnCopy<string> mystr2(mystr1); // reinit
    SimTK_TEST(mystr2 == "unknown" && mystr2.getReinitValue() == "unknown");

    ReinitOnCopy<const char*> mychar1 = "charstr";
    SimTK_TEST(string(mychar1) == "charstr");
    mychar1 = "here is a new string";
    SimTK_TEST(string(mychar1) == "here is a new string");

    ReinitOnCopy<const char*> mychar2(mychar1); // copy construction; reinit
    SimTK_TEST(string(mychar2) == "charstr" 
               && string(mychar2.getReinitValue()) == "charstr");

    ReinitOnCopy<std::vector<char>> vc = {'a','b','c'};
    SimTK_TEST(vc.size()==3 && vc.getReinitValue().size()==3);
    vc = {'x', 'y'};
    SimTK_TEST(vc.size()==2 && vc.getReinitValue().size()==3);

    auto vc2(vc); // copy construction; reinit
    SimTK_TEST(vc2.size()==3 && vc2.getReinitValue().size()==3);
}


int main() {
    SimTK_START_TEST("TestCloneOnWritePtr");
        SimTK_SUBTEST(testEmpty);
        SimTK_SUBTEST(testAllocate);
        SimTK_SUBTEST(testForLeaks);

        SimTK_SUBTEST(testNullOnCopyUniquePtr);
        SimTK_SUBTEST(testResetOnCopy);
        SimTK_SUBTEST(testReinitOnCopy);
    SimTK_END_TEST();
}

