
#ifndef __cdsIomanip_h__
#define __cdsIomanip_h__

#include <sthead.h>

#ifdef USE_STD_STREAMS

#include <iomanip>
using std::setw;
using std::setprecision;
using std::setfill;
using std::ws;
using std::setiosflags;
using std::resetiosflags;

#else

#include <iomanip.h>

#endif

#endif /* __cdsIomanip_h__ */
