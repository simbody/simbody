
#include "cdsRc_ptr.h"

#ifdef TESTING
#include <cdsIomanip.h>
#include "cdsString.h"

#include <assert.h>




template <class T>
int
CDS::rc_ptr<T>::test()
{

// class IntClass {
//   int check;
//   int data;
// public:
//   IntClass(const int& i=0) : data(i), check(1059) {}
//   ~IntClass() { assert(check==1059); }
//   
//   operator int() {return data;}
//
// };

 int exit=0;
 cout << "testing rc_ptr...";
 cout.flush();

 // int size=10;
 rc_ptr<IntClass> a(new IntClass); //default constructor

 *a=20;
 {
   if ( *a != 20 ) {
     cerr << "failure from default constructor" << endl;
     exit=1;
   }
 }
  

 { // this clause to check that destruction works ok.
   rc_ptr<IntClass> b(a);              //copy constructor
 }

 rc_ptr<IntClass> b(a);              //copy constructor
 {
   if ( *b != 20 ) {
     cerr << "failure from copy constructor" << endl;
     exit=1;
   }
 }


 rc_ptr<IntClass> c = b;            //assignment operator
 {
   if ( *c != 20 ) {
     cerr << "failure from assignment operator" << endl;
     exit=1;
   }
 }

 //  not explicitly tested
 // 
 //  T& operator*() const;
 //  operator T*() {return data;};
//
//  T* get() const;
//  T* release();
//  T* reset(T *p = 0);



 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

#endif /* TESTING */
