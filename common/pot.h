#ifndef __pot_hh__
#define __pot_hh__

#include "energyReport.h"
#include <cdsString.h>
#include <cdsRc_ptr.h>
#include <accessor.h>

class DerivList;

/**
 * Virtual base class for potential energy terms.
 */
class Pot {
public: 
    Pot(const String& potName, const String& instanceName)
      : potName_(potName), instanceName_(instanceName), scale_(1.)
    { }
    virtual ~Pot() {}

    //
    // equals operator
    //

    // FIX: is this needed now?
    //bool operator==(const Pot& anotherPot) const {
    //
    //  if ((anotherPot.potName() == potName_) && 
    //      (anotherPot.instanceName() == instanceName_))
    //    return 1;
    //  else
    //    return 0;
    //}

    //
    // accessors
    //

    const char* potName()      const { return potName_; }
    const char* instanceName() const { return instanceName_; }

    //
    // constant, multiplicative scale factor
    //
    float_type scale() const { return scale_; }
    void       setScale(const float_type& x) { scale_ = x; }

    //
    // calculate energy without calculating derivs
    //

    virtual EnergyReport calcEnergy() = 0;

    //
    // calculate energy and derivatives wrt atomic positions. 
    // Adds these derivs to the values in the DerivList that's passed in.
    // 

    virtual EnergyReport calcEnergyAndDerivs(DerivList&) = 0;

    //  Pot* pot() { return this; }

protected:
    //
    // instance vbls
    //

    String potName_;
    String instanceName_;
    float_type scale_;
};

class rc_Pot : public CDS::rc_ptr<Pot> {
public:
    rc_Pot(Pot* ptr=0) : CDS::rc_ptr<Pot>(ptr) {}
};

template<class DerivedPot>
class rc_DerivedPot : public rc_Pot {
public:
    rc_DerivedPot(Pot* realPtr = 0) : rc_Pot(realPtr) {}
    rc_DerivedPot(const rc_DerivedPot& rhs) : rc_Pot(rhs) {}

    //rc_DerivedPot& operator=(const rc_DerivedPot& rhs) { rc_pot = rhs.rc_pot; }
    DerivedPot* operator->() const { return (DerivedPot*)this->ptr(); }
    DerivedPot& operator*() const { return *( (DerivedPot*)this->ptr() ); }

    //operator rc_Pot() { return rc_Pot(*this); }
};

#endif /* __pot_hh__ */
