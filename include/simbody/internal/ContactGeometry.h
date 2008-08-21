#ifndef SimTK_SIMBODY_CONTACT_GEOMETRY_H_
#define SimTK_SIMBODY_CONTACT_GEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/*
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"

namespace SimTK {

class ContactGeometryImpl;

/**
 */
class SimTK_SIMBODY_EXPORT ContactGeometry {
public:
    class HalfSpace;
    class Sphere;
    class HalfSpaceImpl;
    class SphereImpl;
    ContactGeometry() : impl(0) {
    }
    ContactGeometry(const ContactGeometry& src);
    explicit ContactGeometry(ContactGeometryImpl* impl) : impl(impl) {
    }
    virtual ~ContactGeometry();
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;
    ContactGeometry& ContactGeometry::operator=(const ContactGeometry& src);
    bool hasImpl() const {
        return impl != 0;
    }
    const ContactGeometryImpl& getImpl() const {
        assert(impl);
        return *impl;
    }
    ContactGeometryImpl& updImpl() {
        assert(impl);
        return *impl;
    }
    const std::string& getType() const;
    int getTypeIndex() const;
private:
    ContactGeometryImpl* impl;
};

class ContactGeometry::HalfSpace : public ContactGeometry {
public:
    HalfSpace();
};

class ContactGeometry::Sphere : public ContactGeometry {
public:
    Sphere(Real radius);
    Real getRadius() const;
    void setRadius(Real radius);
    const SphereImpl& getImpl() const {
        return static_cast<const SphereImpl&>(getImpl());
    }
    SphereImpl& updImpl() {
        return static_cast<SphereImpl&>(updImpl());
    }
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_GEOMETRY_H_
