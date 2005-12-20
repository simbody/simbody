
#ifndef __cdsExcept_hh__
#define __cdsExcept_hh__

#include <sthead.h>
// mirror of stl exception tree

CDS_BEGIN_NAMESPACE

class exception {
  enum {messLength=300};
public:
  char mess[messLength];

  inline exception(const char* m="");      //constructor
  inline exception(const exception& e);    //copy constructor
  // is assignment operator needed?

#ifdef TESTING
  static int test();
#endif  
};

struct out_of_range : public exception 
{ out_of_range(const char* m="") : exception(m) {} };
struct bad_alloc    : public exception
{ bad_alloc(const char* m="") : exception(m) {} };
//exception
struct SingularError : public exception
{ SingularError(const char* m="SingularError") : exception(m) {} };

inline
exception::exception(const char* m)
{
 int i=0;
 for (const char* p=m ; *p!='\0'&&i<messLength-1 ; p++,i++) 
   mess[i] = *p;
 mess[i] = '\0';
} /* constructor */
 
inline
exception::exception(const exception& e)
{
 int i=0;
 for (const char* p=e.mess ; *p!='\0'&&i<messLength-1 ; p++,i++) 
   mess[i] = *p;
 mess[i] = '\0';
} /* copy constructor */

CDS_END_NAMESPACE


#endif /* __cdsExcept_hh__ */
