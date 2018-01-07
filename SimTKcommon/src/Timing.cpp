/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
 * Defines any routines that have to be supplied to implement the features
 * promised in Timing.h. Currently this consists of faking up the
 * Posix timing routines when using the Microsoft Visual Studio C/C++
 * compiler cl on Windows, or when running on Mac OSX 10.11 or older.
 */


#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Timing.h"

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
#elif SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 
    #include <unistd.h>
    #include <sys/time.h>
    #include <mach/mach.h>
    #include <mach/mach_time.h>
#endif

// These local symbols are not needed on all platforms. Define them just
// when needed to avoid "unused variable" warnings.
#if defined(_MSC_VER) || SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME
    // There are a billion (1e9) nanoseconds in a second.
    static const long long NsPerSec = 1000000000LL;
#endif

#if defined(_MSC_VER)
    // There are a million (1e6) nanoseconds in a millisecond.
    static const long long NsPerMs = 1000000LL;

    // There are a thousand (1e3) nanoseconds in a microsecond.
    static const long long NsPerUs = 1000LL;

    // There are 10 million (1e7) 100ns (hecto ns) ticks in a second.
    static const long long HectoNsPerSec = 10000000LL;

    // There are a million (1e6) microseconds in a second.
    static const long long UsPerSec = 1000000LL;

    // This is the time interval from the Win32 epoch 1/1/1601 to
    // the Unix epoch 1/1/1970, measured in 100ns (hecto ns) ticks.
    // The number of days between those dates (measuring them both
    // using the same calendar) is 134774=(369 years)*365 + 89 leap days,
    // or 11644473600 seconds (=134774*24*3600).
    static const long long FILETIME_1970 = 11644473600LL*HectoNsPerSec;
#endif

// static getnstimeofday()
// -----------------------
// Returns the time of day in seconds and nanoseconds since 1-1-1970.
// This is derived from code by Rolf Steenge, posted at
// https://projects.coin-or.org/Cbc/ticket/45.
// This method is the guts of CLOCK_REALTIME.

#if defined(_MSC_VER)
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
        ts.tv_sec  = (long)  (hecto / HectoNsPerSec);            // seconds
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
        GetSystemTimeAsFileTime(&ft);     // 100-nanoseconds since 1-1-1601
        filetimeToTimespec(ft, *tp);     // now in ns since 1-1-1970
        return 0;
    }
#elif SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 
    static int getnstimeofday(struct timespec *tp) {
        if (!tp) return 0;
        struct timeval tod;
        gettimeofday(&tod, 0); // don't care about timezone
        tp->tv_sec = tod.tv_sec;
        tp->tv_nsec = 1000*tod.tv_usec;
        return 0;
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
#elif SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 
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
#elif SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 
    static int getprocesscputime(struct timespec* tp) {
        if (!tp) return 0;

        task_basic_info        deadThreadTimes; // for terminated threads
        task_thread_times_info liveThreadTimes; // for live threads
        mach_msg_type_number_t deadsz, livesz;
        task_t currentProcess;
        kern_return_t status;
        currentProcess = mach_task_self();
        deadsz = sizeof(deadThreadTimes)/sizeof(int);
        livesz = sizeof(liveThreadTimes)/sizeof(int);
        status = task_info(currentProcess, TASK_BASIC_INFO,
            (task_info_t)&deadThreadTimes,&deadsz);
        status = task_info(currentProcess, TASK_THREAD_TIMES_INFO,
            (task_info_t)&liveThreadTimes,&livesz);

        tp->tv_sec = (time_t)(  deadThreadTimes.user_time.seconds
                              + deadThreadTimes.system_time.seconds
                              + liveThreadTimes.user_time.seconds
                              + liveThreadTimes.system_time.seconds);
        tp->tv_nsec = (long)(1000*(deadThreadTimes.user_time.microseconds
                                   + deadThreadTimes.system_time.microseconds
                                   + liveThreadTimes.user_time.microseconds
                                   + liveThreadTimes.system_time.microseconds));
        while (tp->tv_nsec >= NsPerSec) {
            ++tp->tv_sec;
            tp->tv_nsec -= NsPerSec;
        }
        return 0;
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
#elif SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 
    static int getthreadcputime(struct timespec* tp) {
        if (!tp) return 0;

        thread_basic_info info;
        mach_msg_type_number_t infosz;
        mach_port_t currentThread;
        kern_return_t status;
        currentThread = mach_thread_self(); // get "send rights"
        infosz = sizeof(info)/sizeof(int);
        status = thread_info(currentThread, THREAD_BASIC_INFO,
            (thread_info_t)&info,&infosz);
        // relinquish "send rights"
        mach_port_deallocate(mach_task_self(), currentThread);

        tp->tv_sec  = (time_t)(info.user_time.seconds + info.system_time.seconds);
        tp->tv_nsec = (long)(1000*(  info.user_time.microseconds
                                   + info.system_time.microseconds));
        while (tp->tv_nsec >= NsPerSec) {
            ++tp->tv_sec;
            tp->tv_nsec -= NsPerSec;
        }

        return 0;
    }
#endif

// int clock_gettime()
// -------------------
// Now define the Posix clock_gettime() function in terms of the above helpers.
#if defined(_MSC_VER) || SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 
    int clock_gettime (clockid_t clock_id, struct timespec *tp) {
        int retval = EINVAL;

        switch (clock_id)
        {
        case CLOCK_REALTIME:
            retval = getnstimeofday(tp);
            break;
        case CLOCK_MONOTONIC:
        #ifdef CLOCK_MONOTONIC_HR
        case CLOCK_MONOTONIC_HR:  // "high resolution"; not defined on macOS if
                                  // using SDK MacOSX10.12.sdk or greater with
                                  // deployment target 10.11 or lower.
        #endif
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

