
#ifndef __atomSel_hh__
#define __atomSel_hh__

#include <atom.h>
#include <cdsIostream.h>

class Simulation;
namespace AtomSelAction { class Base; }

//************************************************************
//                   AtomSel class definition
//************************************************************

#include <cdsList.h>
#include <cdsString.h>


class AtomSel {
  Simulation* sim_;
  int ok_;

  String sel;
  CDSList<int> aList_;
  mutable CDSList<bool> boolList_;

  //  AtomSel(const AtomSel&);
public:
  AtomSel(const char*       selStr=""          ,
		Simulation* defaultSimulation=0);

  AtomSel(const CDSList<int>& indices,
		Simulation*   defaultSimulation=0);

  int ok()       const {return ok_;} //FIX: not really supported
  const Simulation* simulation() const {return sim_;}
  Simulation* simulation() {return sim_;}

  const CDSList<int>& indices() const {return aList_;}
  CDSList<int> indices() {return aList_;}
  void setIndices(const CDSList<int>& indices);

  bool containsIndex(const int i) const;
  
  // AtomSel looks like a List of Atoms
  int size() const { return aList_.size(); }
  Atom operator[](const int i) { return Atom(sim_,aList_[i]); }
  const Atom operator[](const int i) const { return Atom(sim_,aList_[i]); }

  String string() const { return sel; }

  //void apply(const AtomSelAction&);
  //void apply(AtomSelAction::Base&);
  void apply(AtomSelAction::Base&) const; //is const a lie?

  friend istream&   operator>>(istream& s, AtomSel& x);
  friend ostream&   operator<<(ostream& s, const AtomSel& x);
};

istream&   operator>>(istream& s, AtomSel& x);
ostream&   operator<<(ostream& s, const AtomSel& x);


#endif /* __atomSel_hh__ */
