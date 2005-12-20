
#ifndef __sthead_hh__
#define __sthead_hh__ 1

#include <cstdlib>

//define FORTRAN macro
#ifndef NO_FORTRAN
#  include <f77defs.h>
#endif /* ! NO_FORTRAN */

/*this is usually nothing at all*/
#define MALLOC_DEBUG

/*this is usually nothing at all*/
#define INSTALL_FPE_HANDLER

/*the old definition*/
#define FMTFLAGS_TYPE long

//architecture-dependent definitions
#include <arch.h>

#ifndef FLOAT_TYPE
#  define FLOAT_TYPE double
#endif
typedef FLOAT_TYPE float_type;


#ifndef TRUE
#     define FALSE 0      
#     define TRUE (!FALSE)
#endif

#ifndef PI
#     define PI 3.14159265358979323846
#endif


#ifndef l_int        /*for the draft standard behavior*/
#    define l_int int
#endif


#ifdef USE_CDS_NAMESPACE
#  define CDS_NAMESPACE(x) CDS::x
#  define CDS_BEGIN_NAMESPACE namespace CDS {
#  define CDS_END_NAMESPACE }
#else
#  define CDS_NAMESPACE(x) x
#  define CDS_BEGIN_NAMESPACE
#  define CDS_END_NAMESPACE
#endif

CDS_BEGIN_NAMESPACE
//really need this?
template<class T> inline
void swap(T &t1,
          T &t2)
{
 T tmp = t1; t1=t2; t2=tmp;
}
CDS_END_NAMESPACE


#ifdef __GNUG__ 
#  if (__GNUC__ > 2)
//stuff for g++-v3
#    define NEWSTYLE_FMTFLAGS
// this flag disables warnings when including the .h files in libstdc++
#    define _CPP_BACKWARD_BACKWARD_WARNING_H
#  endif
#endif

#ifdef NEWSTYLE_FMTFLAGS
#  undef FMTFLAGS_TYPE
#  define FMTFLAGS_TYPE std::ios_base::fmtflags
#endif

#if !__GNUG__
      //really ugly!
#     define istdiostream stdiostream  
#endif

//
// Compaq C++ 6.2 optimizer is a bit too eager: this necessary to 
//  support my use of subVector constructor
//
#ifndef OSF_SUBVECTOR_HACK
#  define OSF_SUBVECTOR_HACK
#endif 


#endif /*__sthead_hh__*/
