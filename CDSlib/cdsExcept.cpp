
#include "cdsExcept.h"


#ifdef TESTING
#include <cdsIostream.h>
#include "cdsString.h"


int
CDS::exception::test()
{
 int exit=1;
 cout << "testing CDSExcept...";
 cout.flush();

 try {
   throw exception("a nice description");
 }
 catch ( CDS::exception e ) {
   if ( String(e.mess) == "a nice description" )
     exit = 0;
   else 
     cout << "mess: >" << e.mess << "< ";
 }

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

#endif /* TESTING */
