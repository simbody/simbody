
#ifndef __dint_xplor_hh__
#define __dint_xplor_hh__

#include "dinternal.h"
#include "xplorVars.h"

class XplorIVM : public IVM {
  int eCount;
  XplorVars xplorVars;
  CDSList<Pair> bondExclude;

public:

  XplorIVM(const XplorVars& xplorVars) : IVM(), eCount(0), xplorVars(xplorVars) {}
  void calcEnergy();

  void syncPos();
  void syncVel();
  void syncDeriv();
  void initXplor();
  ostream& info(ostream&    ostr,
		double&     finalTime,
		const double& stepsize,
		const int&  nprint,
		const int&  nsavc,
		const int&  nsavv,
		const int&  reuseIntVar,
		const int&  steps,
		const int&  nCMresets,
		const char* trajFilename,
		const char* velFilename);
  void outputStepInfo(const int     step,
		      const double& stepsize,
		      const CDSString& type);
  void fixAtomsFromExternal(const int   natom,
			    const long* imove);
  CDSString idAtom(int id) const;
  void entry();
};

#endif /* __dint_xplor_hh__ */
