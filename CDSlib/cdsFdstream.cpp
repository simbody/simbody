#include "cdsFdstream.h"

#ifdef TESTING

#include <cdsString.h>

#include <cdsFstream.h>
#include <stdio.h>

namespace CDS {

int
test_fdstreams()
{
 int exit=0;
 cout << "testing CDS::fdstream...";
 FILE* fp=popen("date","r");

 // open pipe to read from
 if (fp == NULL)
   throw "popen() failed";

 // and initialize input stream to read from it
 fdistream in( fileno(fp) );

 String word;
 in >> word;
 if ( word.length() != 3 ) {
   cerr << "fdstream: error reading first word: " << word << endl;
   exit=1;
 }
 in >> word;
 if ( word.length() != 3 ) {
   cerr << "fdstream: error reading second word: " << word << endl;
   exit=1;
 }

 pclose(fp);

 {
   String tmpFilename = "/tmp/xplor_XXXXXX";
   String testString = "a test string";
   {
     // casting safe here because string length does not change.
     CDS::fdostream oStr( mkstemp((char*)((const char*)tmpFilename)) );
     oStr << testString << '\n';
   }

   FILE* fp=popen("cat " + tmpFilename,"r");
   // open pipe to read from
   if (fp == NULL)
     throw "popen() failed";

   // and initialize input stream to read from it
   fdistream in( fileno(fp) );

   String tmp;
   readline(in,tmp);
   if ( tmp != testString ) {
     cerr << "fdstream: error in file test: " << tmp << endl;
     exit=1;
   }
   pclose(fp);
 }

 cout << (exit?"failed":"ok") << endl;
 return exit;
}

} /* CDS namespace */

#endif /* TESTING */
