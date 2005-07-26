
#ifndef __dihedral_hh__
#define __dihedral_hh__

//
// simple class to evaluate the value of a dihedral angle.
//

#include <atom.h>

class AtomSel;

class Dihedral {
  Atom a0;
  Atom a1;
  Atom a2;
  Atom a3;
  Dihedral(); //inaccessible
public:
  Dihedral(const AtomSel& a0,
	   const AtomSel& a1,
	   const AtomSel& a2,
	   const AtomSel& a3);
  
  float_type value();

  const Atom& atom0() { return a0; }
  const Atom& atom1() { return a1; }
  const Atom& atom2() { return a2; }
  const Atom& atom3() { return a3; }
};

#endif /* __dihedral_hh__ */
