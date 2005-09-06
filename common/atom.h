#ifndef __atom_hh__
#define __atom_hh__

#include "vec3.h"
#include "cdsString.h"

class Simulation;

/**
 * Basic definition for the atom type.  
 *
 * An atom contains information on its name,
 * position/velocity/deriv, and atom-specific
 * parameters like mass, vdw radius, etc.
 */
class Atom {
public:
    static const float_type INVALID_COORD;

    //
    // Constructors
    //

    /* this is only defined because SWIG demands it */
    Atom();

    Atom(      Simulation* sim, int index);
    Atom(const Simulation* sim, int index);

    //
    // accessors
    //

    const char* segmentName() const ;
    const char* residueName() const ;
    const char* atomName() const    ;
    const char* chemType() const    ;
    int   residueNum() const        ;
    Vec3  pos() const;
    Vec3  vel() const;
    const float_type& mass() const  ;
    const float_type& fric() const  ;
    //  const float_type& radius() const;
    const float_type& charge() const;

    void setSegmentName(const char* newName)  ;
    void setResidueName(const char* newName)  ;
    void setAtomName(const char* newName)     ;
    void setChemType(const char* newName)     ;
    void setResidueNum(int newNum)         ;
    void setPos(const Vec3& v);
    void setVel(const Vec3& v);
    void setMass(const float_type& newVal)   ;
    void setFric(const float_type& newVal)   ;
    //  void setRadius(const float_type& newVal) ;
    void setCharge(const float_type& newVal) ;

    String string() const;

    Simulation* simulation() { return sim; }
    const Simulation* simulation() const { return (const Simulation*)sim; }
    int index() const { return index_; }

    // validity test
    bool isValid() const {return (pos().x() < INVALID_COORD);}

protected:
    Simulation* sim;
    int         index_;
};
  
ostream&
operator<<(ostream& s, const Atom& x);

//
// inline function definitions
//


#endif /* __atom_hh__ */
