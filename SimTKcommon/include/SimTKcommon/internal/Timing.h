#ifndef SimTK_SimTKCOMMON_TIMING_H_
#define SimTK_SimTKCOMMON_TIMING_H_

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

/**@file
This file ensures that we have access to the Posix time functions
clock_getttime() and nanosleep(), and also provides some convenient methods
for use in common timing situations. **/

/**@defgroup TimingFunctions Timing Functions
   @ingroup GlobalFunctions

These functions provide a convenient way to do high-precision timing, either 
for real time use or performance measurement, and to sleep for precise time 
intervals. Both elapsed and CPU timing is supported, with the latter on a 
per-process or per-thread basis. Times are measured using a double precision
floating point number of seconds, which is usually the most convenient form.
Alternatives are available that measure time as a 64 bit integer count of 
nanosecond ticks, providing the highest resolution and consistency for very 
short measurements.

Note that you can also use the Posix functions clock_gettime() and nanosleep() 
on any SimTK-supported platform, as well as usleep() (time in microseconds)
although Posix considers that obsolete.

On Linux systems, use of these timing functions requires linking with the 
librt realtime library (-lrt). **/

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"

// This header is needed on Mac and Linux for some or all of the Posix time 
// functions and the timespec struct. We include it on Windows also for 
// uniform cross-platform behavior, since there are many other useful time-
// date-handling symbols declared here on all platforms.
#include <ctime>


// macOS (OSX) 10.12 introduced support for clock_gettime().
// The following logic is based on documentation in /usr/local/Availability.h
// and from this site:
// https://developer.apple.com/library/content/documentation/DeveloperTools/Conceptual/cross_development/Using/using.html.
// This website also recommends using Availability.h over the similarly-helpful
// AvailabilityMacros.h.
#if defined(__APPLE__)
    #include <Availability.h>
#endif
// We'll use the following macros throughout Timing.(h|cpp).
// "MAX_ALLOWED" is the version of the OSX SDK used when building.
// "MIN_REQUIRED" is the "DEPLOYMENT_TARGET": earliest version on which the
//                binaries should run.
// One can use the 10.12 SDK to deploy to earlier releases, like 10.11. The SDK
// version determines if we need to declare clock_gettime() (e.g., developer is
// using an SDK older than 10.12), and the deployment target determines if we
// cannot expect the user's system to contain an implementation of
// clock_gettime() (e.g., user may be running an on OS older than 10.12).
// We need both of these macros because developers may only have the 10.12 SDK,
// but may want to use it to deploy to machines running 10.11. In such a case,
// we cannot declare clock_gettime() ourselves (the SDK does it), but we must
// still define it.
// The number 101200 identifies macOS version 10.12. The explicit version
// number must be used instead of a macro like __MAC_10_12 because pre-10.12
// systems won't have __MAC_10_12 defined.
#if defined(__APPLE__) && __MAC_OS_X_VERSION_MAX_ALLOWED < 101200
    // SDK is older than 10.12.
    #define SimTK_IS_APPLE_AND_MUST_DECLARE_CLOCK_GETTIME 1
#else
    #define SimTK_IS_APPLE_AND_MUST_DECLARE_CLOCK_GETTIME 0
#endif

#if defined(__APPLE__) && __MAC_OS_X_VERSION_MIN_REQUIRED < 101200
    // user's OS may be pre-10.12.
    #define SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 1
#else
    #define SimTK_IS_APPLE_AND_MUST_DEFINE_CLOCK_GETTIME 0
#endif



#if defined(_MSC_VER)
    /* Posix nanosleep() sleeps the indicated number of nanoseconds and returns
    0, or if it is interrupted early it returns how much time was left in 
    rem and returns EINTR. Ours is not interruptable so will always succeed and
    return rem==0. It is OK if rem is NULL, but req==NULL or req<0 returns 
    EINVAL. A time of req==0 is allowed and our interpretation is that the 
    thread relinquishes its time slice to another ready-to-run thread if there
    is one, otherwise returns immediately. This implementation rounds the 
    desired sleep time to the nearest millisecond. On a Linux system, this 
    requires including <time.h> (or <ctime>), which we already included 
    above. */
    SimTK_SimTKCOMMON_EXPORT int nanosleep(const struct timespec* req, struct timespec* rem);

    /* Posix declares this handy function obsolete, but I don't think it is in
    any danger of going away. It sleeps for the given number of microseconds.
    However, using SimTK::sleepInNs() or SimTK::sleepInSec() is safer. */
    typedef unsigned int useconds_t;
    inline int usleep(useconds_t us) {
        struct timespec req;
        req.tv_sec  = (long) (us / 1000000U);
        req.tv_nsec = (long)((us % 1000000U)*1000U);
        int status = nanosleep(&req,0);
        return status ? -1 : 0;
    }
#endif

#if defined(_MSC_VER) || SimTK_IS_APPLE_AND_MUST_DECLARE_CLOCK_GETTIME
    // On Windows and OSX < 10.12, the Posix clock_gettime function is missing.
    typedef long clockid_t;

    /* These constants are the clock ids we support. All the varieties of 
    CLOCK_MONOTONIC are high resolution with no NTP adjustments. I measured 
    the resolutions on a single Windows 7 machine; hopefully they are typical
    (resolution here means how often they are updated):
      - MONOTONIC (counter):    0.001ms      1us
      - REALTIME (time of day):    1ms    1000us
      - CPUTIME (either):         20ms   20000us
    These are slightly conservative resolutions so you should be able to 
    achieve them in practice. */
    #define CLOCK_REALTIME           1 // time of day clock, from 1/1/1970
    #define CLOCK_MONOTONIC          2 // counter from last boot time
    #define CLOCK_MONOTONIC_HR       3 // "high resolution" (same)
    #define CLOCK_MONOTONIC_RAW      4 // "not subject to NTP adjustments" (same)
    #define CLOCK_THREAD_CPUTIME_ID  5 // current thread's cpu time (kernel+user)
    #define CLOCK_PROCESS_CPUTIME_ID 6 // cumulative cpu time of all threads of
                                       //   this process, live or dead

    /* Returns zero if it succeeds (or if tp==NULL); otherwise EINVAL. On a 
    Linux system, this requires including <time.h> (or <ctime>) and linking 
    with -lrt to get the realtime library. */
    SimTK_SimTKCOMMON_EXPORT int clock_gettime(clockid_t clock_id, 
                                               struct timespec *tp); 
#endif



namespace SimTK {

/** @defgroup TimeConversions Timespec/Nanosecond/Second Conversions
    @ingroup TimingFunctions

These inline functions provide a fast and convenient way for doing arithmetic
with the ugly Posix timespec struct. Use them to convert the timespec to
a long long integer number of nanoseconds, do arithmetic in that form, and
then convert back. Negative times are handled correctly (they come up as
the result of subtraction and comparisons). 

We usually prefer to deal with times as a double precision floating point
number of seconds and functions are provided for converting between 
nanoseconds and seconds in this format. 

@par Cautions:
    - For long intervals the precision of time in seconds will necessarily be 
      less than the precision of the nanosecond count, since IEEE double 
      precision has a 53 bit mantissa, while 63 bits are available for the count.
    - A signed long long integer containing a count of nanosecond ticks
      is limited to time intervals of about +/- 292 years, which is 
      substantially less than can be contained in a timespec, but is plenty 
      for interval timing. 
**/

/**@{**/
/** Convert a time stored in a timespec struct to the equivalent number
of nanoseconds (as a signed quantity). @see nsToTimespec() **/    
inline long long timespecToNs(const timespec& ts)
{   return (long long)ts.tv_sec*1000000000LL + (long long)ts.tv_nsec; }

/** Given a signed number of nanoseconds, convert that into seconds and 
leftover nanoseconds in a timespec struct. @see timespecToNs() **/
inline void nsToTimespec(const long long& ns, timespec& ts) 
{   ts.tv_sec  = (long)(ns / 1000000000LL); // signed
    if (ns >= 0) ts.tv_nsec =  (long)(  ns  % 1000000000LL);
    else         ts.tv_nsec = -(long)((-ns) % 1000000000LL); }

/** Given a count of nanosecond ticks as a signed 64 bit integer, return
the same time interval as a double precision floating point number of
seconds. See @ref TimeConversions for cautions. @see secToNs() **/
inline double nsToSec(const long long& ns) 
{   return (double)(ns*SimTK_NS_TO_S); }

/** Given a signed time interval as a double precision floating point number of
seconds, return the same time interval as a count of nanosecond ticks in a 
signed 64 bit integer. See @ref TimeConversions for cautions. @see nsToSec() **/
inline long long secToNs(const double& s) 
{   return (long long)(s*SimTK_S_TO_NS); }
/**@}**/

/** @defgroup CPUTimers Measuring CPU Time
    @ingroup TimingFunctions

These functions provide measurement of CPU time consumed by a process or
by individual threads. Time includes both kernel and user time together,
and is reported as a double precision floating point value in seconds.
CPU timers typically have a very coarse resolution, likely to be in the 
10-50ms range depending on the particulars of your system. That means you 
won't get repeatable results unless you measure substantial amounts of 
CPU time; don't expect to get meaningful information measuring CPU times
of less than a second or so. **/
/**@{**/

/** Return the cumulative CPU time in seconds (both kernel and user time) that 
has been used so far by any of the threads in the currently executing process.

@return CPU time used since this process was created, in seconds, as
a double precision floating point number.
@see threadCpuTime() **/
inline double cpuTime() 
{   timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return (double)(timespecToNs(ts)*SimTK_NS_TO_S); }

/** Return the total CPU time in seconds (both kernel and user time) that 
has been used so far by the currently executing thread.

@return CPU time used since this thread was created, in seconds, as
a double precision floating point number.
@see cpuTime() **/
inline double threadCpuTime() 
{   timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return (double)(timespecToNs(ts)*SimTK_NS_TO_S); }
/**@}**/

/** @defgroup ElapsedTime High-Resolution Elapsed Time Measurement and Sleep
    @ingroup TimingFunctions

These functions provide access to the system's high resolution interval timer,
and highest precision sleep facility for unscheduling a thread until it is 
needed at a specific later time. The interval timer measures elapsed time from
some arbitrary starting point (typically since the system was last booted). It
is expected that this timer provides very precise measurement of short time 
intervals, but cannot be depended upon to measure very long periods without 
drifting. That is, it may not be synchronized to the system time of day.

Generally it is most convenient to measure intervals as a floating point number
of seconds, however this provides a variable amount of precision as the 
absolute value of the timer increases. For maximum precision, you can obtain 
the timer value as an integer number of nanoseconds instead. You can improve
accuracy by subtacting the integer counts first before converting to seconds. 
The actual resolution is system-dependent, but it should be able to accurately
measure elapsed times of 1ms or less, substantially less on some systems. **/
/**@{**/

/** Return current time on the high-resolution interval timer in nanoseconds, 
as a 64-bit integer count. Generally it is more convenient to use realTime() 
which reports the interval time in seconds instead, but the nanosecond count 
is best for maximum accuracy and consistency.

@return Elapsed nanoseconds since some arbitrary time, as a 64 bit integer 
count.
@see realTime(), nsToSec() **/
inline long long realTimeInNs() {
    timespec ts;
    #ifdef CLOCK_MONOTONIC_RAW
        clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    #else
        clock_gettime(CLOCK_MONOTONIC, &ts);
    #endif
    return timespecToNs(ts);
}

/** Return current time on the high-resolution interval timer in seconds. For
maximum precision, you can improve repeatability and accuracy somewhat by 
obtaining the interval times as integer counts using realTimeInNs(). 

@return Elapsed seconds since some arbitrary time, as a double precision 
floating point number.
@see realTimeInNs() **/ 
inline double realTime() {return nsToSec(realTimeInNs());}

/** Sleep for the indicated number of nanoseconds, with the actual precision
system dependent but intended to be the best achievable, hopefully less than
5ms in all cases. However, when this wakes up you should not assume you 
reached the desired time; you could wake up a little earlier or a lot later.
Alternatives to this are %SimTK's sleepInSec() function, or the Posix 
nanosleep() and usleep() functions. 
@see sleepInSec() **/
inline void sleepInNs(const long long& ns) 
{   timespec ts;
    nsToTimespec(ns, ts);
    nanosleep(&ts, 0); }

/** Sleep for the indicated number of seconds, with the actual precision
system dependent but intended to be the best achievable, hopefully less than
5ms in all cases. However, when this wakes up you should not assume you reached
the desired time; you could wake up a little earlier or a lot later. 
Alternatives to this are the %SimTK sleepInNs() function, or the Posix 
nanosleep() and usleep() functions. 
@see sleepInNs() **/
inline void sleepInSec(const double& seconds) {sleepInNs(secToNs(seconds));}
/**@}**/

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_TIMING_H_
