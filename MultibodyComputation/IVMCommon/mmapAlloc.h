
#ifndef __mmapAlloc_hh__
#define __mmapAlloc_hh__

#include <sys/types.h>

class MMapAlloc {
public:

  //
  // the second argument specifies whether to delete the mmaped file.
  //
  //
  // the third argument can be used to allow different processes to attach 
  // to a specified map.
  //

  static void* alloc(size_t,
		     bool deleteFile=1,
		     int  mapID=-1);
  static void free(void*);
};

#endif /* __mmapAlloc_hh__ */
