
#include "cdsAlloc.h"

//
// things to know: 
//
//  new must always return a unique pointer, even when given a size-zero
//  argument
//
//  delete must work when given a NULL-valued pointer (i.e. do nothing)
//
//

//void*
//CDS::DefaultAlloc::alloc(size_t size)
//{
// return new char[size];
//}
//
//void
//CDS::DefaultAlloc::free(void* p)
//{
// delete [] (char*)p;
//}
//
