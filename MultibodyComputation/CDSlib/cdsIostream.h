
#ifndef __cdsIostream_h__
#define __cdsIostream_h__

#include <sthead.h>

#ifdef USE_STD_STREAMS

#define min(a,b) BLAHBLAH_min(a,b)

#include <iostream>
#include <istream>

#undef  min

using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::ends;
using std::ostream;
using std::istream;
using std::iostream;
using std::streambuf;
using std::ios;
using std::flush;

#else

#include <iostream.h>

#endif

#endif /* __cdsIostream_h__ */
