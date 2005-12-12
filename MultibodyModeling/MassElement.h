#ifndef SIMTK_MASS_ELEMENT_H_
#define SIMTK_MASS_ELEMENT_H_

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
 * User-visible definitions for the mass element features which can
 * be owned by bodies.
 */

#include "Feature.h"

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"
#include <iostream>

namespace simtk {

// These are defined here.
class MassElement;
class PointMassElement;
class SphereMassElement;
class CylinderMassElement;
class BrickMassElement;
class GeneralMassElement;

/**
 * This is an abstract interface implemented by any Subsystem
 * that carries mass.
 */
class MassElement : public Feature {
public:
    const RealMeasure&    getMassMeasure() const;
    const StationMeasure& getCentroidMeasure() const;

    static bool               isInstanceOf(const Subsystem&);
    static const MassElement& downcast(const Subsystem&);
    static MassElement&       downcast(Subsystem&);
};

/**
 * A concrete mass element consisting only of a point mass. This has
 * a single parameter "mass" which takes a scalar value, and the point
 * mass element feature requires a station placement on a body.
 */
class PointMassElement : public MassElement {
public:
    explicit PointMassElement(const String&);
    PointMassElement(const PointMassElement&);    // placements are not copied
    PointMassElement& operator=(const PointMassElement&);
    ~PointMassElement();

    // This constructor gives the "mass" parameter a constant value.
    PointMassElement(const String&, const Real&);

    // These create constant placements owned by the PointMassElement
    // feature itself.
    void setMass(const Real&);

    static bool                    isInstanceOf(const Subsystem&);
    static const PointMassElement& downcast(const Subsystem&);
    static PointMassElement&       downcast(Subsystem&);
};

/**
 * A concrete mass element consisting only of a sphere of uniform
 * mass density (not a point!). This has two scalar parameters, "mass"
 * and "radius", and the sphere mass element feature requires a
 * station placement on a body, at which its center will be located.
 */
class SphereMassElement : public MassElement {
public:
    explicit SphereMassElement(const String&);
    SphereMassElement(const SphereMassElement&);    // placements are not copied
    SphereMassElement& operator=(const SphereMassElement&);
    ~SphereMassElement();

    // These create constant placements owned by the SphereMassElement
    // feature itself.
    void setMass    (const Real&);
    void setRadius  (const Real&);

    static bool                    isInstanceOf(const Subsystem&);
    static const PointMassElement& downcast(const Subsystem&);
    static PointMassElement&       downcast(Subsystem&);
};

/**
 * A concrete mass element consisting of a cylinder of uniformly-distributed
 * mass. There are three parameters: mass, radius, and halfLength. The
 * cylinder mass element requires a station and a direction (or plane)
 * placement. The cylinder's center will be located at the station placement,
 * and its long axis aligned with the direction placement. We refer to the
 * end pointed to by the direction as the "top" and the other end as the "bottom"
 * of the cylinder. You can also visualize this as the cylinder's centerline 
 * circle aligned with the plane with the "top" above the plane and "bottom" below.
 */
class CylinderMassElement : public MassElement {
public:
    explicit CylinderMassElement(const String&);
    CylinderMassElement(const CylinderMassElement&);    // placements are not copied
    CylinderMassElement& operator=(const CylinderMassElement&);
    ~CylinderMassElement();

    // These create constant placements owned by the CylinderMassElement
    // feature itself.
    void setMass      (const Real&);
    void setRadius    (const Real&);
    void setHalfLength(const Real&);
    void placeCenter  (const Vec3&);
    void placeAxis    (const Vec3&);

    static bool                       isInstanceOf(const Subsystem&);
    static const CylinderMassElement& downcast(const Subsystem&);
    static CylinderMassElement&       downcast(Subsystem&);
};


} // namespace simtk

#endif // SIMTK_MASS_ELEMENT_H_
