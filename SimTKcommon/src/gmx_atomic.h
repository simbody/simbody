/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
*
* Copyright (c) 2004-2008, Erik Lindahl <lindahl@cbr.su.se>
*
*  Unfortunately, some of the constructs in this file are _very_ sensitive
*  to compiler optimizations and architecture changes. If you find any such
*  errors, please send a message to lindahl@cbr.su.se to help us fix the
*  upstream version too.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*
* And Hey:
* Gnomes, ROck Monsters And Chili Sauce
*/
#ifndef _GMX_ATOMIC_H_
#define _GMX_ATOMIC_H_

/*! \file gmx_atomic.h
 *
 *  @brief Atomic operations for fast SMP synchronization
 *
 *  This file defines atomic integer operations and spinlocks for 
 *  fast synchronization in performance-critical regions of gromacs.
 *
 *  In general, the best option is to use functions without explicit 
 *  locking, e.g. gmx_atomic_fetch_add() or gmx_atomic_cmpxchg().
 *
 *  Not all architecture support atomic operations though inline assembly,
 *  and even if they do it might not be implemented here. In that case
 *  we use a fallback mutex implementation, so you can always count on
 *  the function interfaces working in Gromacs.
 *
 *  Don't use spinlocks in non-performance-critical regions like file I/O.
 *  Since they always spin busy they would waste CPU cycles instead of 
 *  properly yielding to a computation thread while waiting for the disk.
 *
 *  Finally, note that all our spinlock operations are defined to return
 *  0 if initialization or locking completes successfully.
 *  This is the opposite of some other implementations, but the same standard
 *  as used for pthread mutexes. So, if e.g. are trying to lock a spinlock,
 *  you will have gotten the lock if the return value is 0.
 * 
 *  gmx_spinlock_islocked(x) obviously still returns 1 if the lock is locked,
 *  and 0 if it is available, though...
 */



#include <stdio.h>

#include <pthread.h>

#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif




#if ( ( (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__)) && \
        (defined(i386) || defined(__x86_64__)) )                                      \
      || defined (DOXYGEN) )

/* This code is executed for x86 and x86-64, with these compilers:
 * GNU
 * Intel 
 * Pathscale
 * All these support GCC-style inline assembly. 
 * We also use this section for the documentation.
 */

/*! \brief Memory barrier operation
 *
 *  Modern CPUs rely heavily on out-of-order execution, and one common feature
 *  is that load/stores might be reordered. Also, when using inline assembly
 *  the compiler might already have loaded the variable we are changing into
 *  a register, so any update to memory won't be visible.
 *
 *  This command creates a memory barrier, i.e. all memory results before
 *  it in the code should be visible to all memory operations after it - the
 *  CPU cannot propagate load/stores across it.
 */
#define gmx_atomic_memory_barrier() __asm__ __volatile__("": : :"memory")

/* Only gcc and Intel support this check, otherwise set it to true (skip doc) */
#if (!defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined DOXYGEN)
#define __builtin_constant_p(i) (1)
#endif


/*! \brief Gromacs atomic operations datatype
 *
 *  Portable synchronization primitives like mutexes are effective for
 *  many purposes, but usually not very high performance.
 *  One of the problem is that you have the overhead of a function call,
 *  and another is that Mutexes often have extra overhead to make the
 *  scheduling fair. Finally, if performance is important we don't want
 *  to suspend the thread if we cannot lock a mutex, but spin-lock at 100%
 *  CPU usage until the resources is available (e.g. increment a counter).
 *
 *  These things can often be implemented with inline-assembly or other
 *  system-dependent functions, and we provide such functionality for the
 *  most common platforms. For portability we also have a fallback 
 *  implementation using a mutex for locking.
 *
 *  Performance-wise, the fastest solution is always to avoid locking 
 *  completely (obvious, but remember it!). If you cannot do that, the
 *  next best thing is to use atomic operations that e.g. increment a
 *  counter without explicit locking. Spinlocks are useful to lock an
 *  entire region, but leads to more overhead and can be difficult to
 *  debug - it is up to you to make sure that only the thread owning the
 *  lock unlocks it!
 *
 *  You should normally NOT use atomic operations for things like 
 *  I/O threads. These should yield to other threads while waiting for 
 *  the disk instead of spinning at 100% CPU usage.
 *
 *  It is imperative that you use the provided routines for reading
 *  and writing, since some implementations require memory barriers before
 *  the CPU or memory sees an updated result. The structure contents is
 *  only visible here so it can be inlined for performance - it might
 *  change without further notice.
 *
 *  \note No initialization is required for atomic variables.
 *
 *  Currently, we have (real) atomic operations for:
 *
 *  - x86 or x86_64, using GNU compilers
 *  - x86 or x86_64, using Intel compilers 
 *  - x86 or x86_64, using Pathscale compilers
 *  - Itanium, using GNU compilers 
 *  - Itanium, using Intel compilers
 *  - Itanium, using HP compilers
 *  - PowerPC, using GNU compilers 
 *  - PowerPC, using IBM AIX compilers 
 *  - PowerPC, using IBM compilers >=7.0 under Linux or Mac OS X.
 */
typedef struct gmx_atomic
{
    volatile int       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;



/*! \brief Gromacs spinlock
 *
 *  Spinlocks provide a faster synchronization than mutexes,
 *  although they consume CPU-cycles while waiting. They are implemented
 *  with atomic operations and inline assembly whenever possible, and
 *  otherwise we use a fallback implementation where a spinlock is identical
 *  to a mutex (this is one of the reasons why you have to initialize them).
 *
 *  There are no guarantees whatsoever about fair scheduling or
 *  debugging if you make a mistake and unlock a variable somebody
 *  else has locked - performance is the primary goal of spinlocks.
 *
 */
typedef struct gmx_spinlock
{
    volatile unsigned int   lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;





/*! \brief Spinlock static initializer
 *
 *  This is used for static spinlock initialization, and has the same
 *  properties as GMX_THREAD_MUTEX_INITIALIZER has for mutexes.
 *  This is only for inlining in the gmx_thread.h header file. Whether
 *  it is 0, 1, or something else when unlocked depends on the platform.
 *  Don't assume anything about it. It might even be a mutex when using the
 * fallback implementation!
 */
#define GMX_SPINLOCK_INITIALIZER   { 1 }



/*! \brief Return value of an atomic integer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable to read
 *  \return     Integer value of the atomic variable
 */
#define gmx_atomic_read(a)  ((a)->value) 

 
/*! \brief Write value to an atomic integer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable
 *  \param  i   Integer to set the atomic variable to.
 */
#define gmx_atomic_set(a,i)  (((a)->value) = (i))

 
/*! \brief Add integer to atomic variable
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param a   atomic datatype to modify
 *  \param i   integer to increment with. Use i<0 to subtract atomically.
 *
 *  \return The new value (after summation).
 */
static inline int
gmx_atomic_add_return(gmx_atomic_t *     a, 
                      volatile int       i)
{
    int __i;
    
    __i = i;
    __asm__ __volatile__("lock ; xaddl %0, %1;"
                         :"=r"(i) :"m"(a->value), "0"(i));
    return i + __i;
}  
  

/*! \brief Add to variable, return the old value.
 *
 *  This operation is quite useful for synchronization counters.
 *  By performing a fetchadd with N, a thread can e.g. reserve a chunk 
 *  with the next N iterations, and the return value is the index
 *  of the first element to treat.
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param a   atomic datatype to modify
 *  \param i   integer to increment with. Use i<0 to subtract atomically.
 *
 *  \return    The value of the atomic variable before addition.
 */
static inline int
gmx_atomic_fetch_add(gmx_atomic_t *     a,
                     volatile int       i)
{
    int __i;

    __i = i;
    __asm__ __volatile__("lock ; xaddl %0, %1;"
                         :"=r"(i) :"m"(a->value), "0"(i));
    return i;
}


/*! \brief Atomic compare-exchange operation
 *
 *   The \a old value is compared with the memory value in the atomic datatype.
 *   If the are identical, the atomic type is updated to the new value, 
 *   and otherwise left unchanged. 
 *  
 *   This is a very useful synchronization primitive: You can start by reading
 *   a value (without locking anything), perform some calculations, and then
 *   atomically try to update it in memory unless it has changed. If it has
 *   changed you will get an error return code - reread the new value
 *   an repeat the calculations in that case.
 *
 *   \param a        Atomic datatype ('memory' value)
 *   \param oldval   Integer value read from the atomic type at an earlier point
 *   \param newval   New value to write to the atomic type if it currently is
 *                   identical to the old value.
 *
 *   \return The value of the atomic memory variable in memory when this 
 *           instruction was executed. This, if the operation succeeded the
 *           return value was identical to the \a old parameter, and if not
 *           it returns the updated value in memory so you can repeat your
 *           operations on it. 
 *
 *   \note   The exchange occurred if the return value is identical to \a old.
 */
static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *    a, 
                   int               oldval,
                   int               newval)
{
    volatile unsigned long prev;
    
    __asm__ __volatile__("lock ; cmpxchgl %1,%2"
                         : "=a"(prev)
                         : "q"(newval), "m"(a->value), "0"(oldval)
                         : "memory");
    
    return prev;
}


/*! \brief Initialize spinlock
 *
 *  In theory you can call this from multiple threads, but remember
 *  that we don't check for errors. If the first thread proceeded to
 *  lock the spinlock after initialization, the second will happily
 *  overwrite the contents and unlock it without warning you.
 *
 *  \param x      Gromacs spinlock pointer.
 */
static inline void
gmx_spinlock_init(gmx_spinlock_t *   x)
{
    x->lock = 1;
}



/*! \brief Acquire spinlock
 *
 *  This routine blocks until the spinlock is available, and
 *  the locks it again before returning.
 *
 *  \param x     Gromacs spinlock pointer
 */
static inline void
gmx_spinlock_lock(gmx_spinlock_t *  x)
{
    __asm__ __volatile__("\n1:\t" 
                         "lock ; decb %0\n\t" 
                         "jns 3f\n" 
                         "2:\t" 
                         "rep;nop\n\t" 
                         "cmpb $0,%0\n\t" 
                         "jle 2b\n\t" 
                         "jmp 1b\n" 
                         "3:\n\t" 
                         :"=m" (x->lock) : : "memory"); 
}


/*! \brief Attempt to acquire spinlock
 *
 * This routine acquires the spinlock if possible, but if 
 * already locked it return an error code immediately.
 *
 *  \param x     Gromacs spinlock pointer
 *
 * \return 0 if the mutex was available so we could lock it,
 *         otherwise a non-zero integer (1) if the lock is busy.
 */
static inline int
gmx_spinlock_trylock(gmx_spinlock_t *  x)
{
    char old_value;
    
    __asm__ __volatile__("xchgb %b0,%1"
                         :"=q" (old_value), "=m" (x->lock)
                         :"0" (0) : "memory");
    return (old_value <= 0);
}


/*! \brief Release spinlock
 *
 *  \param x     Gromacs spinlock pointer
 *
 *  Unlocks the spinlock, regardless if which thread locked it.
 */
static inline void
gmx_spinlock_unlock(gmx_spinlock_t *  x)
{
    char old_value = 1;
    
    __asm__ __volatile__(
                         "xchgb %b0, %1" 
                         :"=q" (old_value), "=m" (x->lock) 
                         :"0" (old_value) : "memory"
                         );
}
 

/*! \brief Check if spinlock is locked
 *
 *  This routine returns immediately with the lock status.
 *
 *  \param x  Gromacs spinlock pointer
 *
 *  \return 1 if the spinlock is locked, 0 otherwise.
 */
static inline int
gmx_spinlock_islocked(gmx_spinlock_t *  x)
{
    return (*(volatile signed char *)(&(x)->lock) <= 0);
}


/*! \brief Wait for a spinlock to become available
 *
 *  This routine blocks until the spinlock is unlocked, 
 *  but in contrast to gmx_spinlock_lock() it returns without 
 *  trying to lock the spinlock.
 *
 *  \param x  Gromacs spinlock pointer
 */
static inline void
gmx_spinlock_wait(gmx_spinlock_t *   x)
{
    do 
    {
        gmx_atomic_memory_barrier(); 
    } 
    while(gmx_spinlock_islocked(x));
}


#elif ( defined(__GNUC__) && (defined(__powerpc__) || defined(__ppc__)))
/* PowerPC using proper GCC inline assembly. 
 * Recent versions of xlC (>=7.0) _partially_ support this, but since it is
 * not 100% compatible we provide a separate implementation for xlC in
 * the next section.
 */

/* Compiler-dependent stuff: GCC memory barrier */
#define gmx_atomic_memory_barrier() __asm__ __volatile__("": : :"memory")



typedef struct gmx_atomic
{
    volatile int       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;


typedef struct gmx_spinlock
{
    volatile unsigned int   lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;


#define GMX_SPINLOCK_INITIALIZER   { 0 }


#define gmx_atomic_read(a)   ((a)->value) 
#define gmx_atomic_set(a,i)  (((a)->value) = (i))


static inline int
gmx_atomic_add_return(gmx_atomic_t *    a, 
                      int               i)
{
    int t;
    
    __asm__ __volatile__("1:     lwarx   %0,0,%2\n"
                         "\tadd     %0,%1,%0\n"
                         "\tstwcx.  %0,0,%2 \n"
                         "\tbne-    1b"
                         "\tisync\n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value)
                         : "cc" , "memory");
    return t;
}



static inline int
gmx_atomic_fetch_add(gmx_atomic_t *     a,
                     int                i)
{
    int t;
    
    __asm__ __volatile__("\teieio\n"
                         "1:     lwarx   %0,0,%2\n"                         
                         "\tadd     %0,%1,%0\n"
                         "\tstwcx.  %0,0,%2 \n"
                         "\tbne-    1b\n"
                         "\tisync\n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value)
                         : "cc", "memory");
    
    return (t - i);    
}


static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *       a,
                   int                  oldval,
                   int                  newval)
{
    int prev;
    
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\tcmpw    0,%0,%3 \n"
                          "\tbne     2f \n"
                          "\tstwcx.  %4,0,%2 \n"
                          "bne-    1b\n"
                          "\tsync\n"
                          "2:\n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), "m" (a->value)
                          : "cc", "memory");
    
    return prev;
}

static inline void
gmx_spinlock_init(gmx_spinlock_t *x)
{
    x->lock = 0;
}



static inline void
gmx_spinlock_lock(gmx_spinlock_t *  x)
{
    unsigned int tmp;
    
    __asm__ __volatile__("\tb      1f\n"
                         "2:      lwzx    %0,0,%1\n"
                         "\tcmpwi   0,%0,0\n"
                         "\tbne+    2b\n"
                         "1:      lwarx   %0,0,%1\n"
                         "\tcmpwi   0,%0,0\n"
                         "\tbne-    2b\n"
                         "\tstwcx.  %2,0,%1\n"
                         "bne-    2b\n"
                         "\tisync\n"
                         : "=&r"(tmp)
                         : "r"(&x->lock), "r"(1)
                         : "cr0", "memory");
}


static inline int
gmx_spinlock_trylock(gmx_spinlock_t *  x)
{
    unsigned int old, t;
    unsigned int mask = 1;
    volatile unsigned int *p = &x->lock;
    
    __asm__ __volatile__("\teieio\n"
                         "1:      lwarx   %0,0,%4 \n"
                         "\tor      %1,%0,%3 \n"
                         "\tstwcx.  %1,0,%4 \n"
                         "\tbne     1b\n"
                         "\tsync\n"
                         : "=&r" (old), "=&r" (t), "=m" (*p)
                         : "r" (mask), "r" (p), "m" (*p)
                         : "cc", "memory");
    
    return ((old & mask) != 0);    
}


static inline void
gmx_spinlock_unlock(gmx_spinlock_t *  x)
{
    __asm__ __volatile__("\teieio\n": : :"memory");
    x->lock = 0;
}


static inline int
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    return ( x->lock != 0);
}


static inline void
gmx_spinlock_wait(gmx_spinlock_t *x)
{
    do 
    {
        gmx_atomic_memory_barrier(); 
    }
    while(gmx_spinlock_islocked(x));
}



#elif ( (defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM))  && \
        (defined(__powerpc__) || defined(__ppc__)))
/* PowerPC using xlC inline assembly. 
 * Recent versions of xlC (>=7.0) _partially_ support GCC inline assembly
 * if you use the option -qasm=gcc but we have had to hack things a bit, in 
 * particular when it comes to clobbered variables. Since this implementation
 * _could_ be buggy, we have separated it from the known-to-be-working gcc
 * one above.
 */

/* memory barrier - no idea how to create one with xlc! */
#define gmx_atomic_memory_barrier()



typedef struct gmx_atomic
{
    volatile int       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;


typedef struct gmx_spinlock
{
    volatile unsigned int   lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;


#define GMX_SPINLOCK_INITIALIZER   { 0 }


#define gmx_atomic_read(a)   ((a)->value) 
#define gmx_atomic_set(a,i)  (((a)->value) = (i))


static inline int
gmx_atomic_add_return(gmx_atomic_t *    a, 
                      int               i)
{
    int t;
    
    __asm__ __volatile__("1:     lwarx   %0,0,%2 \n"
                         "\t add     %0,%1,%0 \n"
                         "\t stwcx.  %0,0,%2 \n"
                         "\t bne-    1b \n"
                         "\t isync \n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value) );
    return t;
}



static inline int
gmx_atomic_fetch_add(gmx_atomic_t *     a,
                     int                i)
{
    int t;
    
    __asm__ __volatile__("\t eieio\n"
                         "1:     lwarx   %0,0,%2 \n"                         
                         "\t add     %0,%1,%0 \n"
                         "\t stwcx.  %0,0,%2 \n"
                         "\t bne-    1b \n"
                         "\t isync \n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value));
    
    return (t - i);    
}


static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *       a,
                   int                  oldval,
                   int                  newval)
{
    int prev;
    
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\t cmpw    0,%0,%3 \n"
                          "\t bne     2f \n"
                          "\t stwcx.  %4,0,%2 \n"
                          "\t bne-    1b \n"
                          "\t sync \n"
                          "2: \n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), "m" (a->value));
    
    return prev;
}

static inline void
gmx_spinlock_init(gmx_spinlock_t *x)
{
    x->lock = 0;
}



static inline void
gmx_spinlock_lock(gmx_spinlock_t *  x)
{
    unsigned int tmp;
    
    __asm__ __volatile__("\t b      1f \n"
                         "2:      lwzx    %0,0,%1 \n"
                         "\t cmpwi   0,%0,0 \n"
                         "\t bne+    2b \n"
                         "1:      lwarx   %0,0,%1 \n"
                         "\t cmpwi   0,%0,0 \n"
                         "\t bne-    2b \n"
                         "\t stwcx.  %2,0,%1 \n"
                         "\t bne-    2b \n"
                         "\t isync\n"
                         : "=&r"(tmp)
                         : "r"(&x->lock), "r"(1));
}


static inline int
gmx_spinlock_trylock(gmx_spinlock_t *  x)
{
    unsigned int old, t;
    unsigned int mask = 1;
    volatile unsigned int *p = &x->lock;
    
    __asm__ __volatile__("\t eieio\n"
                         "1:      lwarx   %0,0,%4 \n"
                         "\t or      %1,%0,%3 \n"
                         "\t stwcx.  %1,0,%4 \n"
                         "\t bne     1b \n"
                         "\t sync \n"
                         : "=&r" (old), "=&r" (t), "=m" (*p)
                         : "r" (mask), "r" (p), "m" (*p));
    
    return ((old & mask) != 0);    
}


static inline void
gmx_spinlock_unlock(gmx_spinlock_t *  x)
{
    __asm__ __volatile__("\t eieio \n");
    x->lock = 0;
}


static inline void
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    return ( x->lock != 0);
}


static inline void
gmx_spinlock_wait(gmx_spinlock_t *   x)
{
    
    do 
    {
        gmx_atomic_memory_barrier();
    }
    while(gmx_spinlock_islocked(x));
}




#elif (defined(__ia64__) && (defined(__GNUC__) || defined(__INTEL_COMPILER)))
/* ia64 with GCC or Intel compilers. Since we need to define everything through
* cmpxchg and fetchadd on ia64, we merge the different compilers and only provide 
* different implementations for that single function. 
* Documentation? Check the gcc/x86 section.
*/


typedef struct gmx_atomic
{
    volatile int       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;


typedef struct gmx_spinlock
{
    volatile unsigned int   lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;


#define GMX_SPINLOCK_INITIALIZER   { 0 }


#define gmx_atomic_read(a)   ((a)->value) 
#define gmx_atomic_set(a,i)  (((a)->value) = (i))



/* Compiler thingies */
#ifdef __INTEL_COMPILER
void __memory_barrier(void);
int _InterlockedCompareExchange(volatile int *dest, int xchg, int comp);
unsigned __int64 __fetchadd4_rel(unsigned int *addend, const int increment);
/* ia64 memory barrier */
#  define gmx_atomic_memory_barrier() __memory_barrier()
/* ia64 cmpxchg */
#  define gmx_atomic_cmpxchg(a, oldval, newval) _InterlockedCompareExchange(&a->value,newval,oldval)
/* ia64 fetchadd, but it only works with increments +/- 1,4,8,16 */
#  define gmx_ia64_fetchadd(a, inc)  __fetchadd4_rel(a, inc)

#elif defined __GNUC__  
/* ia64 memory barrier */
#  define gmx_atomic_memory_barrier() asm volatile ("":::"memory")
/* ia64 cmpxchg */
static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *   a,
                   int              oldval,
                   int              newval)
{
    volatile int res;
    asm volatile ("mov ar.ccv=%0;;" :: "rO"(oldval));
    asm volatile ("cmpxchg4.acq %0=[%1],%2,ar.ccv":                    
                  "=r"(res) : "r"(&a->value), "r"(newval) : "memory"); 
                          
    return res;
}


/* fetchadd, but on ia64 it only works with increments +/- 1,4,8,16 */
#define gmx_ia64_fetchadd(a, inc)                                             \
({  unsigned long res;                                                        \
    asm volatile ("fetchadd4.rel %0=[%1],%2"                                  \
                  : "=r"(res) : "r"(a), "i" (inc) : "memory");                \
                  res;                                                        \
})


#else /* Unknown compiler */
#  error Unknown ia64 compiler (not GCC or ICC) - modify gmx_thread.h!
#endif



static inline int
gmx_atomic_add_return(gmx_atomic_t *       a, 
                      volatile int         i)
{
    volatile int oldval,newval;    
    volatile int __i = i;

    /* Use fetchadd if, and only if, the increment value can be determined
     * at compile time (otherwise this check is optimized away) and it is
     * a value supported by fetchadd (1,4,8,16,-1,-4,-8,-16).
     */                         
    if (__builtin_constant_p(i) &&
        ( (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) ||         
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) ) )
    {
        oldval = gmx_ia64_fetchadd(a,__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = gmx_atomic_read(a);
            newval = oldval + i;
        }
        while(gmx_atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return newval;
}



static inline int
gmx_atomic_fetch_add(gmx_atomic_t *     a,
                     volatile int       i)
{
    volatile int oldval,newval;    
    volatile int __i = i;
    
    /* Use ia64 fetchadd if, and only if, the increment value can be determined
     * at compile time (otherwise this check is optimized away) and it is
     * a value supported by fetchadd (1,4,8,16,-1,-4,-8,-16).
     */                         
    if (__builtin_constant_p(i) &&
        ( (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) ||         
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) ) )
    {
        oldval = gmx_ia64_fetchadd(a,__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = gmx_atomic_read(a);
            newval = oldval + i;
        }
        while(gmx_atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return oldval;
}


static inline void
gmx_spinlock_init(gmx_spinlock_t *x)
{
    x->lock = 0;
}


static inline void
gmx_spinlock_lock(gmx_spinlock_t *   x)
{
    gmx_atomic_t *a = (gmx_atomic_t *) x;
    unsigned long value;                                                 
    value = gmx_atomic_cmpxchg(a, 0, 1);                             
    if (value)                                                           
    {                                                                    
        do                                                               
        {                                                                
            while (a->value != 0)                                                 
            {                                                            
                gmx_atomic_memory_barrier();                             
            }                                                            
            value = gmx_atomic_cmpxchg(a, 0, 1);                       
        }                                                                
        while (value);                                                   
    }                                                                    
} 


static inline int
gmx_spinlock_trylock(gmx_spinlock_t *   x)
{
    return (gmx_atomic_cmpxchg((gmx_atomic_t *)x, 0, 1) != 0);
}


static inline void
gmx_spinlock_unlock(gmx_spinlock_t *   x)
{
    do
    {
        gmx_atomic_memory_barrier(); 
        x->lock = 0;
    } 
    while (0);
}


static inline int
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    return (x->lock != 0);
}


static inline void
gmx_spinlock_wait(gmx_spinlock_t *   x)
{
    
    do 
    {
        gmx_atomic_memory_barrier();
    }
    while(gmx_spinlock_islocked(x));
}


#undef gmx_ia64_fetchadd



#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
/* HP compiler on ia64 */
#include <machine/sys/inline.h>

#define gmx_atomic_memory_barrier() _Asm_mf()

#define gmx_hpia64_fetchadd(a, i)                           \
    _Asm_fetchadd((_Asm_fasz)_FASZ_W,(_Asm_sem)_SEM_REL,    \
                  (UInt32*)a,(unsigned int) i,              \
                  (_Asm_ldhint)LDHINT_NONE)
 

typedef struct gmx_atomic
{
    volatile int       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;


typedef struct gmx_spinlock
{
    volatile unsigned int   lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;


static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *   a,
                   int              oldval,
                   int              newval)
{
    int ret;
    
    _Asm_mov_to_ar((_Asm_app_reg)_AREG_CCV,(Uint32)oldval,                  
                   (_Asm_fence)(_UP_CALL_FENCE | _UP_SYS_FENCE |         
                                _DOWN_CALL_FENCE | _DOWN_SYS_FENCE));
                   
    ret = _Asm_cmpxchg((_Asm_sz)SZ_W,(_Asm_sem)_SEM_ACQ,(Uint32*)a,    
                       (Uint32)newval,(_Asm_ldhint)_LDHINT_NONE);
                   
    return ret;
}



#define GMX_SPINLOCK_INITIALIZER   { 0 }


#define gmx_atomic_read(a)   ((a)->value) 
#define gmx_atomic_set(a,i)  (((a)->value) = (i))


static inline void 
gmx_atomic_add_return(gmx_atomic_t *       a, 
                      int                  i)
{
    int old,new;    
    int __i = i;
    
    /* On HP-UX we don't know any macro to determine whether the increment
     * is known at compile time, but hopefully the call uses something simple
     * like a constant, and then the optimizer should be able to do the job.
     */                         
    if (  (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) ||         
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) )
    {
        oldval = gmx_hpia64_fetchadd(a,__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = gmx_atomic_read(a);
            newval = oldval + i;
        }
        while(gmx_atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return newval;
}



static inline int
gmx_atomic_fetch_add(gmx_atomic_t *     a,
                     int                i)
{
    int oldval,newval;    
    int __i = i;
    
    /* On HP-UX we don't know any macro to determine whether the increment
     * is known at compile time, but hopefully the call uses something simple
     * like a constant, and then the optimizer should be able to do the job.
     */                         
    if (  (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) ||         
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) )
    {
        oldval = gmx_hpia64_fetchadd(a,__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = gmx_atomic_read(a);
            newval = oldval + i;
        }
        while(gmx_atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return oldval;
}


static inline void
gmx_spinlock_init(gmx_spinlock_t *x)
{
    x->lock = 0;
}





static inline void
gmx_spinlock_trylock(gmx_spinlock_t *x)
{
    int rc;

    rc = _Asm_xchg((_Asm_sz)_SZ_W, (unsigned int *)x, 1        
                    (_Asm_ldhit)_LDHINT_NONE);
    
    return ( (rc>0) ? 1 : 0);
}


static inline void
gmx_spinlock_lock(gmx_spinlock_t *x)
{
    int      status = 1;
    
    do
    {
        if( *((unsigned int *)x->lock) == 0 ) 
        {
            status = gmx_spinlock_trylock(x);
        }
    } while( status != 0);
}


static inline void
gmx_spinlock_unlock(gmx_spinlock_t *   x)
{
    _Asm_fetchadd((_Asm_fasz)_SZ_W,(_Asm_sem)_SEM_REL,                  
                  (unsigned int *)x,-1,(_Asm_ldhint)_LDHINT_NONE);
}



static inline void
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    return ( x->lock != 0 );
}



static inline void
gmx_spinlock_wait(gmx_spinlock_t *   x)
{
    do
    {
        gmx_atomic_memory_barrier(); 
    } 
    while(gmx_spinlock_islocked(x));
}


#undef gmx_hpia64_fetchadd



#elif (defined(_MSC_VER) && (_MSC_VER >= 1200))
/* Microsoft Visual C on x86, define taken from FFTW who got it from Morten Nissov */

#include <windows.h>

#define gmx_atomic_memory_barrier()


typedef struct gmx_atomic
{
    LONG volatile       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;


typedef struct gmx_spinlock
{
    LONG volatile      lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;


#define GMX_SPINLOCK_INITIALIZER   { 0 }




#define gmx_atomic_read(a)  ((a)->value) 
#define gmx_atomic_set(a,i)  (((a)->value) = (i))




#define gmx_atomic_fetch_add(a, i)  \
    InterlockedExchangeAdd((LONG volatile *)a, (LONG) i)

#define gmx_atomic_add_return(a, i)  \
    ( i + InterlockedExchangeAdd((LONG volatile *)a, (LONG) i) )

#define gmx_atomic_cmpxchg(a, oldval, newval) \
    InterlockedCompareExchange((LONG volatile *)a, (LONG) newval, (LONG) oldval)


# define gmx_spinlock_lock(x)   \
    while((InterlockedCompareExchange((LONG volatile *)&x, 1, 0))!=0)


#define gmx_spinlock_trylock(x)   \
    InterlockedCompareExchange((LONG volatile *)&x, 1, 0)


static inline void
gmx_spinlock_unlock(gmx_spinlock_t *   x)
{
    x->lock = 0;
}


static inline int
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    return (*(volatile signed char *)(&(x)->lock) != 0);
}


static inline void
gmx_spinlock_wait(gmx_spinlock_t *   x)
{
    while(gmx_spinlock_islocked(x))
    {
        Sleep(0);
    }
}



#elif defined(__xlC__) && defined (_AIX)
/* IBM xlC compiler on AIX */
#include <sys/atomic_op.h>


#define gmx_atomic_memory_barrier()


typedef struct gmx_atomic
{
    volatile int       value;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_atomic_t;


typedef struct gmx_spinlock
{
    volatile unsigned int      lock;      /*!< Volatile, to avoid compiler aliasing */
}
gmx_spinlock_t;


static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *    a,
                   int               oldval,
                   int               newval)
{
    int t;
    
    if(__check_lock((atomic_p)&a->value, oldval, newval))
    {
        /* Not successful - value had changed in memory. Reload value. */
        t = a->value;
    }
    else
    {
        /* replacement succeeded */
        t = oldval;
    }
    return t;        
}


static inline void 
gmx_atomic_add_return(gmx_atomic_t *       a, 
                      int                  i)
{
    int oldval,newval;    
    
    do
    {
        oldval = gmx_atomic_read(a);
        newval = oldval + i;
    }
    while(__check_lock((atomic_p)&a->value, oldval, newval));

    return newval;
}



static inline void 
gmx_atomic_fetch_add(gmx_atomic_t *       a, 
                     int                  i)
{
    int oldval,newval;    
    
    do
    {
        oldval = gmx_atomic_read(a);
        newval = oldval + i;
    }
    while(__check_lock((atomic_p)&a->value, oldval, newval));
    
    return oldval;
}


static inline void
gmx_spinlock_init(gmx_spinlock_t *   x)
{
    __clear_lock((atomic_p)x,0);
}


static inline void
gmx_spinlock_lock(gmx_spinlock_t *   x)
{
    do
    {
        ;
    }
    while(__check_lock((atomic_p)x, 0, 1));
}


static inline void
gmx_spinlock_trylock(gmx_spinlock_t *   x)
{
    /* Return 0 if we got the lock */
    return (__check_lock((atomic_p)x, 0, 1) != 0)
}


static inline void
gmx_spinlock_unlock(gmx_spinlock_t *   x)
{
    __clear_lock((atomic_p)x,0);
}


static inline void
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    return (*((atomic_p)x) != 0);
}


static inline void
gmx_spinlock_wait(gmx_spinlock_t *    x)
{
    while(gmx_spinlock_islocked(x)) { ; } 
}


#else
/* No atomic operations, use mutex fallback. Documentation is in x86 section */


#define gmx_atomic_memory_barrier()

/* System mutex used for locking to guarantee atomicity */
static pthread_mutex_t
gmx_atomic_mutex = PTHREAD_MUTEX_INITIALIZER;


typedef struct gmx_atomic
{
    int       value;
}
gmx_atomic_t;

#define gmx_spinlock_t     pthread_mutex_t

 
#  define GMX_SPINLOCK_INITIALIZER   PTHREAD_MUTEX_INITIALIZER

/* Since mutexes guarantee memory barriers this works fine */
#define gmx_atomic_read(a)   ((a)->value)


static inline void
gmx_atomic_set(gmx_atomic_t *   a, 
               int              i)
{
    /* Mutexes here are necessary to guarantee memory visibility */
    pthread_mutex_lock(&gmx_atomic_mutex);
    a->value = i;
    pthread_mutex_unlock(&gmx_atomic_mutex);
}


static inline int
gmx_atomic_add_return(gmx_atomic_t *   a, 
                      int              i)
{
    int t;
    pthread_mutex_lock(&gmx_atomic_mutex);
    t = a->value + i;
    a->value = t;
    pthread_mutex_unlock(&gmx_atomic_mutex);
    return t;
}


static inline int
gmx_atomic_fetch_add(gmx_atomic_t *   a,
                     int              i)
{
    int old_value;
    
    pthread_mutex_lock(&gmx_atomic_mutex);
    old_value  = a->value;
    a->value   = old_value + i;
    pthread_mutex_unlock(&gmx_atomic_mutex);
    return old_value;
}


static inline int
gmx_atomic_cmpxchg(gmx_atomic_t *           a, 
                   int                      oldv,
                   int                      newv)
{
    int t;
    
    pthread_mutex_lock(&gmx_atomic_mutex);
    t = a->value;
    if (t == oldv)
    {
        a->value = newv;
    }
    pthread_mutex_unlock(&gmx_atomic_mutex);
    return t;
}


#define gmx_spinlock_init(lock)       pthread_mutex_init(lock)
#define gmx_spinlock_lock(lock)       pthread_mutex_lock(lock)
#define gmx_spinlock_trylock(lock)    pthread_mutex_trylock(lock)
#define gmx_spinlock_unlock(lock)     pthread_mutex_unlock(lock)

static inline int
gmx_spinlock_islocked(gmx_spinlock_t *   x)
{
    int rc;
    
    if(gmx_spinlock_trylock(x) != 0)
    {
        /* It was locked */
        return 1;
    }
    else
    {
        /* We just locked it */
        gmx_spinlock_unlock(x);
        return 0;
    }
}


static inline void
gmx_spinlock_wait(gmx_spinlock_t *   x)
{
    int rc;
    
    gmx_spinlock_lock(x);
    /* Got the lock now, so the waiting is over */
    gmx_spinlock_unlock(x);
}


#endif




/*! \brief Spinlock-based barrier type
 *
 *  This barrier has the same functionality as the standard
 *  gmx_thread_barrier_t, but since it is based on spinlocks
 *  it provides faster synchronization at the cost of busy-waiting.
 *
 *  Variables of this type should be initialized by calling
 *  gmx_spinlock_barrier_init() to set the number of threads
 *  that should be synchronized.
 */
typedef struct gmx_spinlock_barrier
{
    gmx_atomic_t            count;     /*!< Number of threads remaining     */
    int                     threshold; /*!< Total number of threads         */
    volatile int            cycle;     /*!< Current cycle (alternating 0/1) */
}
gmx_spinlock_barrier_t;
 



/*! \brief Initialize spinlock-based barrier
 *
 *  \param barrier  Pointer to _spinlock_ barrier. Note that this is not
 *                  the same datatype as the full, thread based, barrier.
 *  \param count    Number of threads to synchronize. All threads
 *                  will be released after \a count calls to 
 *                  gmx_spinlock_barrier_wait().  
 */
static inline void 
gmx_spinlock_barrier_init(gmx_spinlock_barrier_t *         barrier,
                          int                              count)
{
    barrier->threshold = count;
    barrier->cycle     = 0;
    gmx_atomic_set(&(barrier->count),count);
}




/*! \brief Perform busy-waiting barrier synchronization
*
*  This routine blocks until it has been called N times,
*  where N is the count value the barrier was initialized with.
*  After N total calls all threads return. The barrier automatically
*  cycles, and thus requires another N calls to unblock another time.
*
*  Note that spinlock-based barriers are completely different from
*  standard ones (using mutexes and condition variables), only the 
*  functionality and names are similar.
*
*  \param barrier  Pointer to previously create barrier.
*
*  \return The last thread returns -1, all the others 0.
*/
static inline int
gmx_spinlock_barrier_wait(gmx_spinlock_barrier_t *   barrier)
{
  int    cycle;
  int    status;
  
  /* We don't need to lock or use atomic ops here, since the cycle index 
    * cannot change until after the last thread has performed the check
    * further down. Further, they cannot reach this point in the next 
    * barrier iteration until all of them have been released, and that 
    * happens after the cycle value has been updated.
    *
    * No synchronization == fast synchronization.
    */
  cycle = barrier->cycle;
  
  /* Decrement the count atomically and check if it is zero.
    * This will only be true for the last thread calling us.
    */
  if( gmx_atomic_add_return( &(barrier->count), -1 ) == 0)
  { 
    gmx_atomic_set(&(barrier->count), barrier->threshold);
    barrier->cycle = !barrier->cycle;
    
    status = -1;
  }
  else
  {
    /* Wait until the last thread changes the cycle index.
    * We are both using a memory barrier, and explicit
    * volatile pointer cast to make sure the compiler
    * doesn't try to be smart and cache the contents.
    */
    do
    { 
      gmx_atomic_memory_barrier();
    } 
    while( *(volatile int *)(&(barrier->cycle)) == cycle);
    
    status = 0;
  }
  return status;
}




#ifdef __cplusplus
}
#endif


#endif /* _GMX_ATOMIC_H_ */
