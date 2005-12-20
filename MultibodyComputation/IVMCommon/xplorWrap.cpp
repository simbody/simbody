
/*
 * convenience C++ wrappers for xplor routines
 */

#ifdef NOTDEF

#include "xplorWrap.h"
#include <xplorSimulation.h>
#include <f77defs.h>
#include "prototypes.h"
#include <cdsExcept.h>
#include <cdsAuto_arr.h>
#include <xplorVars.h>
#include <cdsFdstream.h>

#include <cdsFstream.h>
#include <string.h>
#include <unistd.h>



int
XplorWrap::openFile(const char* fname)
{
 long ret;
 long error;
 long oldWrnlev = xplorVars()->wrnlev;
 if ( xplorVars()->wrnlev <= 5 )  // to get rid of assfil's open message
   xplorVars()->wrnlev = 4;
 FORTRAN(assfil)(fname,ret,"READ","FORMATTED",error,
		 strlen(fname),strlen("READ"),strlen("FORMATTED"));
 xplorVars()->wrnlev = oldWrnlev;
 if ( error )
   throw CDS_NAMESPACE(exception)( String("XplorWrap::openFile: error opening")
				   + " " +fname);
 FORTRAN(instrm)(ret);
 return ret;
} /* openFile */

bool
XplorWrap::closeFile(int unit)
  // close current file
{
 long error;
 FORTRAN(vclose)(unit,"KEEP",error,strlen("KEEP"));
 FORTRAN(destrm)();
 return 1;
} /* closeFile */

CDSList< String > 
XplorWrap::command(const char*            cmd,
		   const CDSList<String>& names)
{
 String tmpFilename = "/tmp/xplor_XXXXXX";
 {
   // casting safe here because string length does not change.
   CDS::fdostream oStr( mkstemp((char*)((const char*)tmpFilename)) );
   oStr << cmd << '\n';
   oStr << "end\n" << ends;
   //oStr.close(); //shouldn't be necessary
   //              // under alpha/linux, sometimes file is not closed
   //              // before the following openFile...
 }
 int unit = openFile(tmpFilename);

 simulation->syncTo(); // only necessary fields

 FORTRAN(xplor_parse)(1);

 simulation->syncFrom(); 
 closeFile(unit);
 unlink(tmpFilename);

 CDSList< String > ret;
 for (int i=0 ; i<names.size() ; i++) {
   String var = names[i];
   var.upcase();
   if (var[0]!='$') 
     var = "$" + var;
   char val[1024];
   long vallen;
   FORTRAN(getvar_asstring)(var,val,vallen,1024,var.length(),sizeof(val));
   ret.append( String(val,vallen) );
 }
 return ret;
} /* command */

CDSList< int > 
XplorWrap::select(const char* selString)
  //
  // note that xplor selections are offset 1
  //
{
 simulation->syncTo(); // in case we changed something

 String tmpFilename = "/tmp/xplor_XXXXXX";
 {
   // casting safe here because string length does not change.
   CDS::fdostream oStr( mkstemp((char*)((const char*)tmpFilename)) );
   oStr << '(' << selString << ')' << '\n' << ends;
   //oStr.close(); //shouldn't be necessary
 }
 int unit = openFile(tmpFilename);

 long oldEcho   = xplorVars()->echo;
 long oldWrnlev = xplorVars()->wrnlev;
 xplorVars()->echo = 0;
 if ( xplorVars()->wrnlev <= 5 )  // to get rid of assfil's open message
   xplorVars()->wrnlev = 4;

 auto_arr<long > flags(new long[xplorVars()->natom]);
 long nflags;
 FORTRAN(selcta)(flags,nflags,
		 xplorVars()->x,xplorVars()->y,xplorVars()->z,1);
 FORTRAN(makind)(flags,xplorVars()->natom,nflags); 

 xplorVars()->wrnlev = oldWrnlev;
 xplorVars()->echo   = oldEcho;
 closeFile(unit);
 unlink(tmpFilename);
 
 CDSList<int> ret(nflags,1);
 for (int i=0 ; i<nflags ; i++)
   ret[i] = flags[i]-1;  //correct index offset

 return ret;
} /* select */

void
XplorWrap::shell()
{
 simulation->syncTo(); // only necessary fields

 FORTRAN(xplor_parse)(1);

 simulation->syncFrom();
} /* shell */

void
XplorWrap::setRandomSeed(const double& seed) const //const??
{
 xplorVars()->randomSeed = seed;
}

 
double 
XplorWrap::randomSeed() const
{
 return xplorVars()->randomSeed;
}

double
XplorWrap::uniformRandom() const
{
 double ret=0;
 FORTRAN(ggubfs)(ret);
 return ret;
}

XplorVars* XplorWrap::xplorVars_=0;

XplorVars*
XplorWrap::xplorVars()
{
 return xplorVars_;
}

void
XplorWrap::resetXplorVars(XplorVars* vars)
{
 delete xplorVars_;
 xplorVars_ = vars;
}


extern "C" void
FORTRAN(xplorvars_reset)(const int&          natom,
			       double* const x,
			       double* const y,
			       double* const z,
			       double* const xv,
			       double* const yv,
			       double* const zv,
			       double* const dx,
			       double* const dy,
			       double* const dz,
			       char* const   segid, 
			       char* const   res, 
			       char* const   resid, 
			       char* const   type,
			       char* const   iac,
			 const int&          nenert,
			       long*         qener,    
			 const char* const   aner,    
			 const double*       renr,    
			 const int&          nbond,
			 const long*   	   ib,
			 const long*   	   jb,
			       double*       amass,
			       double*       fbeta,
			       double*       charge,
			       char* const   wd,
			 const long*         imove,
			 const double&       timeFac,
			 const double&       kBoltzman,
			       double&       seed,
			       long&         wrnlev,
			       long&         echo)
{
 XplorWrap::resetXplorVars( new XplorVars(natom,x,y,z,xv,yv,zv,dx,dy,dz,
					  segid, res, resid, type, iac,
					  nenert,qener,aner,renr,nbond,ib,jb,
					  amass,fbeta,charge,imove,wd,
					  timeFac,kBoltzman,seed,wrnlev,echo) );
} /* xplorvars_reset */

#endif
