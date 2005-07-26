
#ifndef __cdsAlloc_hh__
#define __cdsAlloc_hh__

#include <sthead.h>

CDS_BEGIN_NAMESPACE

class DefaultAlloc {
public:

  static void* alloc(size_t size) { return new char[size]; }
  static void free(void* p) { delete [] (char*)p; }
};

CDS_END_NAMESPACE

#endif /* __cdsAlloc_hh__ */
