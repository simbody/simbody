
#ifndef __prototypes_hh__
#define __prototypes_hh__

//
// x-plor subroutines which we call
//

#include <f77defs.h>

extern "C" {
void FORTRAN(xplor_init,XPLOR_INIT)(long&  wrnlen,
			 long&  argc,
			 char** argv);
void FORTRAN(xplor_cleanup,XPLOR_CLEANUP)();
void FORTRAN(xplor_parse,XPLOR_PARSE)(const long& qsubshell);
void FORTRAN(energy,ENERGY)(const double* x,
		     const double* y,
		     const double* z);
void FORTRAN(instrm,INSTRM)(const long& unit);
void FORTRAN(destrm,DESTRM)();
void FORTRAN(assfil,ASSFIL)(const char* name,
				      long&  coordUnit,
				const char* mode,
				const char* form,
				      long&  err,
				const int   lname,
				const int   lmode,
				const int   lform);
void FORTRAN(writtc,WRITTC)(const char* str,
				const long& natom,
				const double* x,
				const double* y,
				const double* z,
				const long& crdnum,
				const long* crdind,
				const long& stepnum,
			              long& qfirsc,
				const double& timeStep,
				const long& nsavc,
				const long& coordUnit,
				const long &formattedOutput,
				const char* format,
				const double& scale,
				const double& offset,
				const int     lstr,
				const int     lformat);
void FORTRAN(declar,DECLAR)(const char*   name,
				const char*   type,
				const char*   stinit,
				const char*   dont_use_me,
				const double& dpinit,
				const int     lname,
				const int     ltype,
				const int     lstinit);
void FORTRAN(pusend,PUSEND)(const char* prompt,
				const int   len);
void FORTRAN(nextwd,NEXTWD)(const char* prompt,
				const int   len);
void FORTRAN(nextst,NEXTST)(const char* prompt,
				      char* ret,
				const int   lprompt,
				const int   lret);
void FORTRAN(miscom,MISCOM)(const char *prompt,
				      long&  used,
				const int   len);
void FORTRAN(nextlo,NEXTLO)(const char *prompt,
				      long&  used,
				const int   len);
void FORTRAN(nexti,NEXTI)(const char *prompt,
				     long&  used,
			       const int   len);
void FORTRAN(nextf,NEXTF)(const char*   prompt,
				     double& used,
			       const int     len);
void FORTRAN(nextfi,NEXTFI)(const char* prompt,
				      char* ret,
				const int lprompt,
				const int lret);
void FORTRAN(chkend,CHKEND)(const char *prompt,
				      long&  used,
				const int   len);
void FORTRAN(selcta,SELCTA)(      long* flags,
				      long& nflags,
				const double* x,
				const double* y,
				const double* z,
				const long& blah);
void FORTRAN(makind,MAKIND)(      long* flags,
				const long& natom,
				const long& nflags); 
void FORTRAN(vclose,VCLOSE)(const long&  unit,
				const char* str,
				      long&  error,
				const int   lstr);
void FORTRAN(printe,PRINTE)(const double* renrl);
void FORTRAN(tdynkin,TDYNKIN)(const long& dof,
				 double* xv,
				 double* yv,
				 double* zv,
			    const double* amass,
			    const long* imove);
void FORTRAN(getvar_asstring,GETVAR_ASSTRING)(const char* var,
				    char* val,
				    long&  vallen,
			      const long&  maxvallen,
			      const int   lvar,
			      const int   lval);
void FORTRAN(ggubfs,GGUBFS)(double& randomNum);
}

#endif /* __prototypes_hh__ */
