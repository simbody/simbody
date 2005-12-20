#include "cdsSStream.h"
#include "cdsMath.h"


#include <cdsIomanip.h>

int
CDSStringStream_test()
{
 int exit=0;
 cout << "testing CDSSStream...";

 {
   StringStream str;
   str << "hello";
   str << ends;
   
   if (str.str() != "hello") {
     cerr << "CDSStringStream: failure: " << str.str() << " != hello" << endl;
     exit=1;
   }
 }

 {
 }

 {
   IStringStream str( "42 hello" );
   int num;
   String s;
   str >> num >> s;
   if ( num != 42 || s != "hello" ) {
     cerr << "CDSIStringSream: failure: " << num << ' ' << s << endl;
     exit=1;
   }
 }

 {
   OStringStream str;
   str << "hello2";
   const char* c = str.str();
   str << " + more!";

   if (str.str() != "hello2 + more!") {
     cerr << "CDSOStringStream: str() failure: >" << str.str() 
	  << "< != >hello2 + more!<" 
	  << endl;
     exit=1;
   }
 }

 {
   OStringStream str;
   str << "hello2";
   const char* c = str.str_volatile();
   str << " + more!";

   if (str.str_volatile() != "hello2 + more!") {
     cerr << "CDSOStringStream: str_volatile() failure: >" 
	  << str.str_volatile() 
	  << "< != >hello2 + more!<" 
	  << endl;
     exit=1;
   }
 }

#ifdef DEBUG_ALLOC
 { // stress tests
   for (int i=0 ; i<20 ; i++) {
     OStringStream s;
     for (int j=0 ; j<100000 ; j++)
       s << "a";
     s << ends;
     for (int j=0 ; j<100000 ; j++)
       if ( s.str()[j] != 'a' ) {
	 cerr << "CDSOStringStream failure in stress test.\n";
	 exit=1;
       }
     
     cout << i << ' ';
     cout.flush();
   }
 }
#endif /* DEBUG_ALLOC */

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

