/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 * Defines any routines that have to be supplied to implement the features
 * promised in Timing.h. Currently this consists of faking up the
 * Posix timing routines when using the Microsoft Visual Studio C/C++ 
 * compiler cl on Windows, or when running on Mac OSX.
 */


#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Timing.h"

#include "pthread.h" // included in SimTK_SDK/include

#include <string>
#include <cstring>
#include <cctype>
#include <ctime>

#include "errno.h"

#if defined(_MSC_VER)
    // Keeps MS VC++ quiet about sprintf, strcpy, etc.
    #pragma warning(disable:4996)

    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
    #include <Windows.h>
#elif defined(__APPLE__)
    #include <unistd.h>
    #include <mach/mach_time.h>
#endif

// There are a billion (1e9) nanoseconds in a second.
static const long long NsPerSec = 1000000000LL;

// There are a million (1e6) nanoseconds in a millisecond.
static const long long NsPerMs = 1000000LL;

// There are a thousand (1e3) nanoseconds in a microsecond.
static const long long NsPerUs = 1000LL;

// There are 10 million (1e7) 100ns (hecto ns) ticks in a second.
static const long long HectoNsPerSec = 10000000LL;

// There are a million (1e6) microseconds in a second.
static const long long UsPerSec = 1000000LL;

// static getnstimeofday()
// -----------------------
// Returns the time of day in seconds and nanoseconds since 1-1-1970.
// This is derived from code by Rolf Steenge, posted at 
// https://projects.coin-or.org/Cbc/ticket/45.
// This method is the guts of CLOCK_REALTIME.
#if defined(_MSC_VER)
    // This is the time interval from the Win32 epoch 1/1/1601 to
    // the Unix epoch 1/1/1970, measured in 100ns (hecto ns) ticks.
    // The number of days between those dates (measuring them both
    // using the same calendar) is 134774=(369 years)*365 + 89 leap days,
    // or 11644473600 seconds (=134774*24*3600).
    static const long long FILETIME_1970 = 11644473600LL*HectoNsPerSec;

    // Return the value of a FILETIME reinterpreted as a
    // long long integer count of the number of 100ns ("hecto ns")
    // ticks since start of 1/1/1601 UTC.
    static inline long long
    filetimeToHectoNs(const FILETIME& ft) {
        // Must do this by copying, not casting, to avoid alignment
        // problems; FILETIME doesn't have to be 16 byte aligned.
        LARGE_INTEGER large;
        large.HighPart = ft.dwHighDateTime;
        large.LowPart = ft.dwLowDateTime;
        return large.QuadPart;
    }

    // Input is number of 100ns ticks since some base time; output is that
    // same interval changed to seconds and nanoseconds in a timespec struct.
    static inline void
    hectoNsToTimespec(const long long& hecto, struct timespec& ts) {
        ts.tv_sec  = (long)  (hecto / HectoNsPerSec);	        // seconds
        ts.tv_nsec = (long) ((hecto % HectoNsPerSec) * 100LL); // nanoseconds
    }

    // Given time reported as a Win32 FILETIME measured from 1/1/1601, convert
    // it to a timespec measured from 1/1/1970.
    // CAUTION: note shift of epoch (measurement base) here.
    static inline void filetimeToTimespec(const FILETIME& ft, 
                                          struct timespec& ts) {
        long long hecto = filetimeToHectoNs(ft);
        hectoNsToTimespec(hecto - FILETIME_1970, ts);
    }

    static int getnstimeofday(struct timespec *tp) {
        if (!tp) return 0;
        FILETIME ft;
        GetSystemTimeAsFileTime(&ft);	 // 100-nanoseconds since 1-1-1601
        filetimeToTimespec(ft, *tp);     // now in ns since 1-1-1970
        return 0;
    }
#elif defined(__APPLE__)
    static int getnstimeofday(struct timespec *tp) {
        if (!tp) return 0;
        //TODO
        tp->tv_sec = 0; tp->tv_nsec = 0;
        return EINVAL;
    }
#endif

// static getperformancecounter()
// ------------------------------
// This is the guts of CLOCK_MONOTONIC (all variants). This is the number of 
// seconds and nanoseconds since an arbitrary start time, with a very high 
// resolution. This is for measuring short intervals very accurately; long 
// term it might drift away from the real time clock.
#if defined(_MSC_VER)
    static int getperformancecounter(struct timespec *tp) {
        if (!tp) return 0;

        // The frequency is not allowed to change after the system is started
        // (according to Microsoft) so we just call this once.
        static long long ticksPerSec = 0LL;
        if (ticksPerSec==0LL) {
            LARGE_INTEGER tps;
            QueryPerformanceFrequency(&tps);
            ticksPerSec = tps.QuadPart;
        }

        LARGE_INTEGER ticks;
        if(!QueryPerformanceCounter(&ticks)) 
            return EINVAL; // bad news

        // Careful: don't try to simplify or rearrange these expressions; 
        // we're depending on integer arithmetic done in long long precision.
        tp->tv_sec = (long) (ticks.QuadPart / ticksPerSec);
        tp->tv_nsec =(long)((ticks.QuadPart % ticksPerSec)*NsPerSec / ticksPerSec);

        return 0;
    }
#elif defined(__APPLE__)
    static int getperformancecounter(struct timespec *tp) {
        if (!tp) return 0;
  
        static mach_timebase_info_data_t info = {0,0};   
  
        if (info.denom == 0)   
                mach_timebase_info(&info);   

        unsigned long long ticks = mach_absolute_time();
        unsigned long long ticksInNs = (ticks * info.numer) / info.denom;   
  
        tp->tv_sec  = (time_t)(ticksInNs / NsPerSec);   
        tp->tv_nsec = (long)(ticksInNs % NsPerSec);  

        return 0;
    }
#endif

// static getprocesscputime()
// --------------------------
// This is the guts of CLOCK_PROCESS_CPUTIME_ID. It returns the number of
// seconds and nanoseconds of cpu time used by all the threads of the
// currently executing process since it started.
#if defined(_MSC_VER)
    static int getprocesscputime(struct timespec* tp) {
        if (!tp) return 0;

        FILETIME creationTime, exitTime, kernelTime, userTime; 
        if (!GetProcessTimes(GetCurrentProcess(),
            &creationTime, &exitTime,
            &kernelTime, &userTime))
            return EINVAL;

        long long ktime = filetimeToHectoNs(kernelTime);
        long long utime = filetimeToHectoNs(userTime);
        hectoNsToTimespec(ktime+utime, *tp);
        return 0;
    }
#elif defined(__APPLE__)
    static int getprocesscputime(struct timespec* tp) {
        if (!tp) return 0;
        //TODO
        tp->tv_sec = 0; tp->tv_nsec = 0;
        return EINVAL;
    }
#endif

// static getthreadcputime()
// -------------------------
// This is the guts of CLOCK_THREAD_CPUTIME_ID. It returns the number of
// seconds and nanoseconds of cpu time used by the currently executing 
// thread since it started.
#if defined(_MSC_VER)
    static int getthreadcputime(struct timespec* tp) {
        if (!tp) return 0;

        FILETIME creationTime, exitTime, kernelTime, userTime; 
        if (!GetThreadTimes(GetCurrentThread(),
            &creationTime, &exitTime,
            &kernelTime, &userTime))
            return EINVAL;

        long long ktime = filetimeToHectoNs(kernelTime);
        long long utime = filetimeToHectoNs(userTime);
        hectoNsToTimespec(ktime+utime, *tp);
        return 0;
    }
#elif defined(__APPLE__)
    static int getthreadcputime(struct timespec* tp) {
        if (!tp) return 0;

        //TODO
        tp->tv_sec = 0; tp->tv_nsec = 0;
        return EINVAL;
    }
#endif

// Now define the Posix clock_gettime() function in terms of the above helpers.
#if defined(_MSC_VER) || defined(__APPLE__)
    int clock_gettime (clockid_t clock_id, struct timespec *tp) {
        int retval = EINVAL;

        switch (clock_id)
        {
        case CLOCK_REALTIME:
            retval = getnstimeofday(tp);
            break;
        case CLOCK_MONOTONIC:
        case CLOCK_MONOTONIC_HR:  // "high resolution"
        case CLOCK_MONOTONIC_RAW: // "not subject to NTP adjustments"
            retval = getperformancecounter(tp);
            break;
        case CLOCK_PROCESS_CPUTIME_ID:
            retval = getprocesscputime(tp);
            break;
        case CLOCK_THREAD_CPUTIME_ID:
            retval = getthreadcputime(tp);
            break;
        default:
            if (tp) {tp->tv_sec = 0; tp->tv_nsec = 0;}
            retval = EINVAL;
            break;
        }

        return retval;
    }
#endif

// int nanosleep()
// ---------------
// Sleep for the amount of time required in req. If we return
// early, rem is supposed to say how much time is left. However, our
// implementation is noninterruptable so we'll always report a
// zero remainder. Note that 0 is an allowable value and just means
// "give up your time slice" which will allow another ready-to-run
// thread to start if there is one but otherwise return immediately.
#if defined(_MSC_VER)
    int nanosleep(const struct timespec* req, struct timespec* rem) {
        if (!req) return EINVAL;
        const long long reqns = SimTK::timespecToNs(*req);
        if (reqns < 0LL) return EINVAL;
        // Round to the nearest ms.
        static const long long halfMsInNs = NsPerMs/2;
        const DWORD reqms = (DWORD)((reqns + halfMsInNs)/NsPerMs);
        Sleep(reqms);
        if (rem) {rem->tv_sec=0; rem->tv_nsec=0;}
        return 0;
    }
#endif

