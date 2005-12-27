
#ifndef __semaphore_hh__
#define __semaphore_hh__

#include <sthead.h>

CDS_BEGIN_NAMESPACE

class Semaphore {
  void* sem;
  //Semaphore(); //default constructor- inaccessible
public:
  //blocked argument:
  //  the initial semaphore value.
  //  0 - if want first call to wait to block
  //  1 - if want first call to succeed, subsequent to block (for wake)
  //      this is apropriate for a mutex
  Semaphore(bool blocked=0);
  ~Semaphore();

  void wait(); //wait until wake
  void wake(); //wake processes in wait()
  int  query(); //return semaphore value
};


CDS_END_NAMESPACE

#endif /* __semaphore_hh__ */
