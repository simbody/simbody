//
// wrapper around xplor commands for scripting language
//

#ifndef __xplorwrap_hh__
#define __xplorwrap_hh__

#include <cdsList.h>
#include <cdsString.h>

class XplorVars;
class Simulation;
class XplorSimulation;

class XplorWrap {
  static XplorVars*  xplorVars_;
  XplorSimulation* simulation;
  int openFile(const char* fname);
  bool closeFile(int unit);

  XplorWrap(); //inaccessible

public:
  XplorWrap(Simulation* simulation) : 
    simulation( (XplorSimulation*)simulation ) {}

  // evaluates an xplor expression
  // names is a list of xplor variable names 
  // Returns a list of the values of those variables in the 
  // same order as they were requested.
  //
  CDSList< String > command(const char*            cmd,
			    const CDSList<String>& names);
  void shell();

  // allows you to eval an xplor selection from python
  //
  CDSList< int > select(const char* selString);

  //random number interface
  void   setRandomSeed(const double& newVal) const;
  double randomSeed() const;
  double uniformRandom() const;

  static XplorVars* xplorVars();
  static void resetXplorVars(XplorVars*);
};

#endif /* __xplorwrap_hh__ */
