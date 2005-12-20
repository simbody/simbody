
#ifndef __xplorVars_hh__
#define __xplorVars_hh__

class XplorVars {
  //  XplorVars& operator=(const XplorVars&);
public:
  //XplorVars(const XplorVars&); //default copy constructor ok

  const int &natom;
  double* const x;
  double* const y;
  double* const z;
  double* const xv;
  double* const yv;
  double* const zv;
  double* const dx;
  double* const dy;
  double* const dz;
  char* const segName;
  char* const resName;
  char* const resNum;    // xplor stores resid as characters
  char* const atomName;
  char* const chemType;
  const int    nenert;
  long*  qener;
  const char* const aner;
  const double* renr;    //energy array
  const int& nbond;  //bond list info
  const long* ib;
  const long* jb;
  double *mass;
  double *fbeta;
  double *charge;
  const long* imove;
  char * const wd;
  const double timeFac;
  const double kBoltzman;
  double& randomSeed;
  long& wrnlev;
  long& echo;
  
  XplorVars(const int&          natom,
		  double *const x,
		  double *const y,
		  double *const z,
	     	  double *xv,
	     	  double *yv,
	     	  double *zv,
		  double *const dx,
		  double *const dy,
		  double *const dz,
                  char*   const segName,
                  char*   const resName,
	          char*   const resNum,
                  char*   const atomName,
                  char*   const chemType,
	    const int&          nenert,
		  long          qener[],
	    const char* const   aner,
	    const double        renr[],
	    const int&          nbond,
	    const long          ib[],
	    const long          jb[],
		  double *mass,
		  double *fbeta,
		  double *charge,
	    const long          imove[],
		  char * const wd,
	    const double& timeFac,
	    const double& kBoltzman,
	          double& seed,
		  long&   wrnlev,
		  long&   echo) :
    natom(natom), x(x), y(y), z(z), xv(xv), yv(yv), zv(zv),
    dx(dx), dy(dy), dz(dz), 
    segName(segName), resName(resName), resNum(resNum), atomName(atomName),
    chemType(chemType),
    nenert(nenert), qener(qener), aner(aner), renr(renr), 
    nbond(nbond), ib(ib),  jb(jb), 
    mass(mass), fbeta(fbeta), charge(charge),
    imove(imove), wd(wd), timeFac(timeFac), 
    kBoltzman(kBoltzman), randomSeed(seed), wrnlev(wrnlev), echo(echo) {}
  
};

//extern XplorVars* xplorVars;


#endif /* __xplorVars_hh__ */
