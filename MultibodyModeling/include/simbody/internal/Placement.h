#ifndef SIMTK_PLACEMENT_H_
#define SIMTK_PLACEMENT_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * Declaration of the client side Placement class, and associated global operators.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Geometry.h"
#include "simbody/internal/Mechanics.h"

#include <iostream>

namespace simtk {

class Subsystem;
class Feature;
class PlacementValue;

/**
 * A Placement is an expression to be used to determine a specification for
 * a particular Feature. In general it involves operators on the placements of
 * other Features. The implementation is opaque.
 */
class Placement {
public:
    Placement() : rep(0) { }
    Placement(const Placement&);
    Placement& operator=(const Placement&);
    ~Placement();

    // Implicit conversions to a Placement of appropriate concrete type.
    Placement(const Real&);
    Placement(const Vec3&);
    Placement(const Mat33&);
    Placement(const UnitVec3&);
    Placement(const RotationMat&);
    Placement(const MatInertia&);
    Placement(const Frame&);

    // These create a "feature placement" (reference to the placement
    // of a feature) of a particular type, but starting with a generic
    // Subsystem, which must turn out to be a Feature. Note that
    // it is not necessary for the Feature to have a placement at the time
    // we create this reference; it need only have one by the time we
    // require a numerical value.
    // The first constructor below is also an implicit
    // conversion, so any Feature may be used as a Placement. The second
    // still references the indicated Feature (not one of its subfeatures)
    // but then selects a subcomponent of the resulting Placement as 
    // indicated by the index.
    Placement(const Subsystem&);
    Placement(const Subsystem&, int index);

    void realize(Stage) const;

    // This will throw an exception if any of this Placement's dependencies
    // have not yet been realized. Note that this will always succeed if
    // this is a constant Placement.
    PlacementValue   calcValue() const;

    bool isConstant()         const; // a plain old numerical value?
    bool isFeatureReference() const; // just a reference to a feature?
    const Feature& getReferencedFeature() const; // only if isFeatureReference()==true

    bool dependsOn(const Feature&) const; // recursive dependency check

    String toString(const String& linePrefix="") const;
    String getPlacementTypeName() const; // for messages only

    // For internal use only.
    explicit Placement(class PlacementRep*);
    bool                      hasRep() const            {return rep != 0;}
    const class PlacementRep& getRep() const            {assert(rep); return *rep;}
    class PlacementRep&       updRep()                  {assert(rep); return *rep;}
    void                      setRep(PlacementRep* pp)  {assert(!rep); rep=pp;}
    void checkPlacementConsistency(const Subsystem* expOwner, 
                                   int              expIndexInOwner,
                                   const Subsystem& expRoot) const;
protected:
    class PlacementRep* rep;
    friend class PlacementRep;
};

// Global operators involving Placements. Note that these actually
// represent families of operators overloaded based on their argument types.

std::ostream& operator<<(std::ostream& o, const Placement&);

// unary
Placement negate    (const Placement&);
Placement abs       (const Placement&);
Placement sqrt      (const Placement&);
Placement square    (const Placement&);
Placement exp       (const Placement&);
Placement log       (const Placement&);
Placement sin       (const Placement&);
Placement cos       (const Placement&);
Placement asin      (const Placement&);
Placement acos      (const Placement&);
Placement length    (const Placement&);
Placement normalize (const Placement&);

// binary
Placement add       (const Placement& l, const Placement& r); 
Placement subtract  (const Placement& l, const Placement& r); 
Placement multiply  (const Placement& l, const Placement& r); 
Placement divide    (const Placement& l, const Placement& r);
Placement distance  (const Placement& l, const Placement& r);
Placement angle     (const Placement& l, const Placement& r);
Placement dot       (const Placement& l, const Placement& r);
Placement cross     (const Placement& l, const Placement& r);

// alternative access to some of above routines via overloaded operators
inline Placement operator- (const Placement& p)                     {return negate(p);}
inline Placement operator+ (const Placement& l, const Placement& r) {return add(l,r);} 
inline Placement operator- (const Placement& l, const Placement& r) {return subtract(l,r);}
inline Placement operator* (const Placement& l, const Placement& r) {return multiply(l,r);}
inline Placement operator/ (const Placement& l, const Placement& r) {return divide(l,r);}

} // namespace simtk

#endif // SIMTK_PLACEMENT_H_
