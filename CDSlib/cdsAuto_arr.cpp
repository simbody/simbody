
#include "cdsAuto_arr.h"

#include "cdsIomanip.h"
#include "cdsString.h"

int
auto_arr_test()
{
 int exit=0;
 cout << "testing auto_arr...";
 cout.flush();

 int size=10;
 auto_arr<int> a(new int[size]); //default constructor

 for (int i=0 ; i<size ; i++)
   a[i] = size-i;
 {
   bool ok=1;
   for (int i=0 ; i<size; i++)
     if ( a[i] != size-i )
       ok=0;
   if ( !ok ) {
     cerr << "failure from default constructor" << endl;
     exit=1;
   }
 }
  

 auto_arr<int> b(a);              //copy constructor
 {
   bool ok=1;
   for (int i=0 ; i<size; i++)
     if ( b[i] != size-i ) ok=0;
   b[3] = 42;
   if ( !ok || b[3] != 42 ) {
     cerr << "failure from copy constructor" << endl;
     exit=1;
   }
 }


 auto_arr<int> c = b;            //assignment operator
 {
   bool ok=1;
   c[3] = 7;
   for (int i=0 ; i<size; i++)
     if ( c[i] != size-i ) ok=0;
   if ( !ok ) {
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
