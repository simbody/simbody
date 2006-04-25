/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
#include <algorithm>
#include <functional>

class MyAbstractType {
public:
    virtual ~MyAbstractType() { }
    virtual SimTK::String myName() const = 0;
    virtual MyAbstractType* clone() const = 0;

    virtual bool operator==(const MyAbstractType&) const = 0;
};

class ConcreteOne : public MyAbstractType {
public:
    ConcreteOne(const SimTK::String s) : objectName(s) { }
    SimTK::String myName() const { return "ConcreteOne::" + objectName; }
    MyAbstractType* clone() const { return new ConcreteOne(*this); }
    bool operator==(const MyAbstractType& a) const {
        return dynamic_cast<const ConcreteOne*>(&a) &&
               objectName == dynamic_cast<const ConcreteOne&>(a).objectName;
    }
private:
    SimTK::String objectName;
}; 

class Num : public MyAbstractType {
public:
    Num(SimTK::Real d) : objectNum(d) { }
    SimTK::String myName() const { return "Num::" + SimTK::String(objectNum); }
    MyAbstractType* clone() const { return new Num(*this); }
    bool operator==(const MyAbstractType& a) const {
        return dynamic_cast<const Num*>(&a) &&
               objectNum == dynamic_cast<const Num&>(a).objectNum;
    }
private:
    SimTK::Real objectNum;
};

std::ostream&
operator<<(std::ostream& o, const MyAbstractType& v)
    { o << v.myName(); return o; }
    
template <class T> std::ostream&
operator<<(std::ostream& o, const std::vector<T>& v)
{ for (size_t i=0; i<v.size(); ++i) {o<<v[i]<<" ";} return o;}

int main()
{
  try {
    SimTK_DEBUG("Running ListTest ...\n");

    int major,minor,build;
    SimTK_version_SimTKcommon(&major,&minor,&build);
    std::printf("SimTKcommon library version: %d.%d.%d\n", major, minor, build);

    char out[100];
    const char* keylist[] = { "version", "library", "type", "debug", "garbage", "authors", "copyright", "svn_revision", 0 };
    std::printf("SimTK_about_SimTKcommon():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKcommon(*p, 100, out);
        std::printf("  about(%s)='%s'\n", *p, out);
    }

    std::printf("PATH=%s\n", getenv("PATH"));

    //SimTK_ASSERT3(1 < 0, "help!!! %d %g %s", 99, 1.23, "just a test");
    //int i=9;
    //try{SimTK_INDEXCHECK(10,i,15,"ouch");} 
    //catch (const std::exception& e) {cout<<e.what()<<endl; }
    // SimTK_SIZECHECK(i+3,8,"justAtest");

    const int ints[] = {1,2,9,-5,2,1};

    std::vector<int> vint;
    vint.assign(ints, ints+6);
    cout << "std::vector: " << vint << endl;
    vint.erase(vint.begin()+2);
    cout << vint << endl;
    std::sort(vint.begin(),vint.end(),std::greater<int>());
    cout << vint << endl;

    SimTK::Array<int> aint(6, ints);
    cout << "SimTK::Array: " << aint << endl;

    std::sort(aint.begin(),aint.end()-3,std::less<int>());
    cout << "aint=" << aint << endl;


    SimTK::ArrayView<int> av = aint(3,3);
    cout << "av=view aint(3,3):" << (SimTK::Array<int>)av << endl;
    cout << "av.front()=" << av.front() << " back=" << av.back() << endl;

    std::sort(av.begin(),av.end());
    cout << "after sort, av=" << (SimTK::Array<int>)av << endl;
    cout << "  and original aint=" << aint << endl;
    
    try {
        
        
    const MyAbstractType* a[] = { new ConcreteOne("hi"), new Num(23.), new Num(-27.2), new ConcreteOne("bye") };
    cout << '|'; for (int i=0; i<4; ++i) cout << a[i] << '|'; cout << endl;

    SimTK::List<MyAbstractType> myl(4,a);
    cout<<"myl{"; for (int i=0; i<myl.size(); ++i) std::cout << '<' << myl[i] << '>'; cout<<'}'<<endl;
    cout << '|'; for (int i=0; i<4; ++i) std::cout << a[i] << '|'; cout << endl;

    // keep Purify happy
    for (int i=0; i<4; ++i) delete a[i], a[i]=0;

    cout << "myl.find(23): " << myl.find(Num(23.)) << endl;

    myl.reverse();
    cout<<"REVERSED: myl{"; for (int i=0; i<myl.size(); ++i) cout << '<' << myl[i] << '>'; cout<<'}'<<endl;
    cout << '|'; for (int i=0; i<4; ++i) cout << &myl[i] << '|'; cout << endl;

    myl[3] = ConcreteOne("this is the replacement");
    myl.push_back(Num(99.29));
    myl[0] = new Num(12344);
    cout<<"myl{"; for (int i=0; i<myl.size(); ++i) cout << '<' << myl[i] << '>'; cout<<'}'<<endl;

    // auto conversion of ListView to List works fine on gcc, like this:
    //   const SimTK::List<double>& dup = myl(3,2);
    // but Microsoft thinks the conversion allows it to destruct the ListView
    SimTK::ListView<MyAbstractType> dup=myl(2,3); 
    SimTK::List<MyAbstractType> copy(dup);   
    std::cout<<"\ncopy{"; for (int i=0; i<copy.size(); ++i) std::cout << '<' << copy[i] << '>'; std::cout<<'}'<<std::endl;
    std::cout<<"myl{"; for (int i=0; i<myl.size(); ++i) std::cout << '<' << myl[i] << '>'; std::cout<<'}'<<std::endl;
    copy[0]=copy[1]=copy[2]=Num(0.001);
    std::cout<<"\ncopy{"; for (int i=0; i<copy.size(); ++i) std::cout << '<' << copy[i] << '>'; std::cout<<'}'<<std::endl;
    std::cout<<"myl{"; for (int i=0; i<myl.size(); ++i) std::cout << '<' << myl[i] << '>'; std::cout<<'}'<<std::endl;

    std::cout<<"\ndup{"; for (int i=0; i<dup.size(); ++i) std::cout << '<' << dup[i] << '>'; std::cout<<'}'<<std::endl;
    dup[0]=new Num(10101.); dup[1]=new Num(20202.); dup[2]=new ConcreteOne("should be a string now!");
    std::cout<<"dup{"; for (int i=0; i<dup.size(); ++i) std::cout << '<' << dup[i] << '>'; std::cout<<'}'<<std::endl;
    std::cout<<"myl{"; for (int i=0; i<myl.size(); ++i) std::cout << '<' << myl[i] << '>'; std::cout<<'}'<<std::endl;
    
    const SimTK::List<MyAbstractType> shl(*(const SimTK::List<MyAbstractType>*)&myl,3,2);
    std::cout<<"\nshl{"; for (int i=0; i<shl.size(); ++i) std::cout << '<' << shl[i] << '>'; std::cout<<'}'<<std::endl;
    //shl[0] = Num(1); shl[1] = new Num(2);
    std::cout<<"shl{"; for (int i=0; i<shl.size(); ++i) std::cout << '<' << shl[i] << '>'; std::cout<<'}'<<std::endl;
    std::cout<<"myl{"; for (int i=0; i<myl.size(); ++i) std::cout << '<' << myl[i] << '>'; std::cout<<'}'<<std::endl;
    
    }
    catch(const SimTK::Exception::Base& b)
        { std::cout << b.getMessage() << std::endl; } 
        
         
    try {
    
    const double ddd[] = { 11, 12, 13, 14, 15, 16 };

    SimTK::Array<double> ad(6, ddd);
    //for (int i=0; i<6; ++i) ad.push_back(ddd[i]);
    
    SimTK_DEBUG("Whole ...\n");
    for (int i=0; i<ad.size(); ++i)
        std::cout << ad[i] << std::endl; 
        
    const SimTK::Array<double> adconst(ad(2,3));
    SimTK_DEBUG("Adconst ...\n");
    for (int i=0; i<adconst.size(); ++i)
        std::cout << adconst[i] << std::endl;   

    // Careful: if we make a non-const view of a const Array, we'll be 
    // allowed to call non-const methods, but the dangerous ones should
    // catch the problem at run time.
    SimTK::Array<double> adnonconst31(adconst,1,2);
    try { 
        std::cout << "adnonconst31[0]=" << adnonconst31[0] << std::endl;
    } catch (const SimTK::Exception::Base& b) {
        std::cout << "Expected exception: Should not allow operator[]:" << std::endl;
        std::cout << b.getMessage() << std::endl;
    } 
    try { 
        adnonconst31 = ad(3,2);
    } catch (const SimTK::Exception::Base& b) {
        std::cout << "Expected exception: Should not allow operator=:" << std::endl;
        std::cout << b.getMessage() << std::endl;
    } 
    
    const SimTK::Array<double> adconst31(adconst,1,2);
    try { 
        std::cout << "adconst31[0]=" << adconst31[0] << std::endl;
    } catch (const SimTK::Exception::Base& b) {
        std::cout << "UNEXPECTED EXCEPTION. Should have allowed operator[]:" << std::endl;
        std::cout << b.getMessage() << std::endl;
    }    
     
    // This shouldn't compile:              
    // adconst31 = ad(3,1);
      
    SimTK_DEBUG("Adconst ...\n");
    for (int i=0; i<adconst.size(); ++i)
        std::cout << adconst[i] << std::endl; 
                
    SimTK_DEBUG("Slice 3,2 of Whole ...\n");
    SimTK::Array<double> adx(ad(3,2));
    for (int i=0; i<adx.size(); ++i)
        std::cout << adx[i] << std::endl; 
    adx[0] = 9999.;
    adx[1] = 12345.;
    SimTK_DEBUG("Slice after assignment...\n");
    for (int i=0; i<adx.size(); ++i)
        std::cout << adx[i] << std::endl; 
    SimTK_DEBUG("Whole(3,2) now ...\n");
    for (int i=0; i<ad(3,2).size(); ++i)
        std::cout << ad(3,2)[i] << std::endl; 
        
{        
    SimTK_DEBUG("VIEW 3,2 of Whole ...\n");
    // auto conversion of ArrayView to Array works fine on gcc, like this:
    //   const SimTK::Array<double>& ddup = ad(3,2);
    // but Microsoft thinks the conversion allows it to destruct the ArrayView
    const SimTK::ArrayView<double>& ddup = ad(3,2);
    const SimTK::ArrayView<double> adx(ddup);
    for (int i=0; i<adx.size(); ++i)
        std::cout << adx[i] << std::endl; 
    ad[3] = 29999.;
    ad[4] = 212345.;
    SimTK_DEBUG("Whole(3,2) after assignment ...\n");
    for (int i=0; i<ad(3,2).size(); ++i)
        std::cout << ad(3,2)[i] << std::endl;
    SimTK_DEBUG("VIEW now...\n");
    for (int i=0; i<adx.size(); ++i)
        std::cout << adx[i] << std::endl; 
         
}

    std::cout << "ELEMENTWISE ASSIGNMENT:" << std::endl;
    SimTK::Array<double> tryInit(10);
    for (int i=0; i<tryInit.size(); ++i)
        std::cout << tryInit[i] << std::endl; 
    tryInit = 1.23;
    for (int i=0; i<tryInit.size(); ++i)
        std::cout << tryInit[i] << std::endl;        
    }
    catch(const SimTK::Exception::Base& b)
    {
        std::cout << b.getMessage() << std::endl;
    }           

    std::cout << "SimTK_LAPACK_INTERFACE=" <<
#if (SimTK_LAPACK_INTERFACE == SimTK)
    "SimTK"
#elif (SimTK_LAPACK_INTERFACE == ACML)
    "ACML"
#else
    "UNKNOWN"
#endif
    << std::endl;

  } catch(const std::exception& e) {
      std::cout << "std::exception: " << e.what() << std::endl;
  }
    return 0;
}


