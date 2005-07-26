#ifndef __f77defs_h__
#define __f77defs_h__

/*
 * all ia32 compilers for linux seem to use this naming convention
 */

#ifdef WIN32
    #define STDCALL __stdcall
    #define FORTRAN(x,X) X
    #define FORTRAN_STRLEN_FOLLOWS_DECL ,int
    #define FORTRAN_STRLEN_FOLLOWS_CALL(n) ,n
    #define FORTRAN_STRLEN_ATEND_DECL
    #define FORTRAN_STRLEN_ATEND_CALL(n)
#else
    #define STDCALL
    #define FORTRAN(x,X) x ## _
    #define FORTRAN_STRLEN_FOLLOWS_DECL
    #define FORTRAN_STRLEN_FOLLOWS_CALL(n)
    #define FORTRAN_STRLEN_ATEND_DECL ,int
    #define FORTRAN_STRLEN_ATEND_CALL(n) ,n
#endif

//#define FORTRAN(x) (x ## _)

#include <g77defs.h>

#endif /* __f77defs_h__ */
