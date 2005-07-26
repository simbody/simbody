
#ifndef __mutex_hh__
#define __mutex_hh__

#include <sthead.h>

CDS_BEGIN_NAMESPACE

class Mutex {
  void *mx; //void* so that it can be used with different schemes
public:
  Mutex();
  ~Mutex();
  void hold();
  void release();
  class Error {};
};

class MutexHold {
  MutexHold(); //default constructor
  MutexHold(const MutexHold&); //copy constructor
  MutexHold &operator=(const MutexHold&); //assignment operator
  
  Mutex &m;
public:
  MutexHold(CDS_NAMESPACE(Mutex)& m) : m(m) { m.hold(); }
  ~MutexHold() { m.release(); }

};

CDS_END_NAMESPACE

#ifdef USE_CDS_NAMESPACE
using CDS::MutexHold;
#endif /* USE_CDS_NAMESPACE */

#endif /* __mutex_hh__ */
