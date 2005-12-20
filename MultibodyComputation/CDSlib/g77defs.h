#ifndef __g77defs_h__
#define __g77defs_h__

#ifdef USING_G77
#undef FORTRAN_INIT

//
// call these routines to intialize the GNU Fortran RTL
//
//

#ifdef CPLUSPLUS
extern "C" {
#endif
  void f_setsig();
  void f_init();
#ifdef CPLUSPLUS
}
#endif

#define FORTRAN_INIT f_setsig(); f_init();
#endif /* USING_G77 */

#endif /* __g77defs_h__ */
