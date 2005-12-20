#ifndef __dint_atom_hh__
#define __dint_atom_hh__

#include "vec3.h"

//#include "atom.h"

class AtomClusterNode;

class IVMAtom {
public:    
    IVMAtom(const int    ix,
            const double &m)
      : index(ix), mass(m), fric(0.), bonds(0,1),
        node(0), station_B(-1.2345e67), station_G(-1.2345e67),
        pos(0.), vel(0.), force(0.)
    {}
/*
    IVMAtom(const int   ix,
            const Atom& atom)
      : index(ix), mass(0.), fric(0.), bonds(0,1),
        node(0), station_B(-1.2345e67), station_G(-1.2345e67),
        pos(0.), vel(0.), force(0.)
    { 
        *this = atom; // sets pos,vel,mass,fric
    }

    IVMAtom& operator=(const Atom& a) { 
        pos = a.pos(); vel = a.vel(); mass = a.mass(); fric = a.fric(); 
        return *this;
    }
*/
public:
    // atom information
    int    index;   // slot # in the atoms list
    double mass;
    double fric;
    
    // molecule information
    CDSList<IVMAtom*>   bonds;

    // clustering information
    AtomClusterNode*    node;
    Vec3                station_B;    // location in node's local frame

    // Initial and then calculated spatial (Cartesian) info.
    Vec3 station_G; // station expressed in G (still measured from OB)
    Vec3 pos;
    Vec3 vel;
    Vec3 force; // this is the Cartesian force applied to this atom
};

ostream& operator<<(ostream &s,const IVMAtom* a);

#endif /* __dint_atom_hh__ */
