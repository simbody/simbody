
#ifndef __cdsFstream_h__
#define __cdsFstream_h__

#include <sthead.h>

#ifdef USE_STD_STREAMS

#include <fstream>

using std::ofstream;
using std::ifstream;

#else

#include <fstream.h>

#endif

#endif /* __cdsFstream_h__ */
