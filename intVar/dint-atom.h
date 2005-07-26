#ifndef __dint_atom_hh__
#define __dint_atom_hh__

#include <atom.h>

class HingeNode;

class IVMAtom {
public:
    int                 index;
    HingeNode*          node;
    CDSList<IVMAtom*>   bonds;

    // Cartesian info.
    Vec3 pos;
    Vec3 vel;
    Vec3 deriv; // this is Cartesian force (i.e., energy gradient), not accel

    double mass;
    double fric;

    IVMAtom(const int index,
            const double &mass)
      : index(index), node(0), bonds(0,1),
        pos(0.), vel(0.), deriv(0.), mass(mass), fric(0.)
    {}

    IVMAtom(const int index,
            const Atom& atom)
      : index(index), node(0), bonds(0,1), deriv(0.)
    { 
        *this = atom; 
    }

    IVMAtom&
    operator=(const Atom& a) { 
        pos = a.pos(); vel = a.vel(); mass = a.mass(); fric = a.fric(); 
        return *this;
    }
};

ostream& operator<<(ostream &s,const IVMAtom* a);

#endif /* __dint_atom_hh__ */
