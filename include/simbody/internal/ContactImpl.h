#ifndef SimTK_SIMBODY_CONTACT_IMPL_H_
#define SimTK_SIMBODY_CONTACT_IMPL_H_

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

#include "simbody/internal/Contact.h"

namespace SimTK {

/**
 * This is the internal implementation class for Contact.
 */
    
class SimTK_SIMBODY_EXPORT ContactImpl {
public:
    ContactImpl(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth);
    virtual ~ContactImpl() {
        assert(referenceCount == 0);
    }
protected:
    friend class Contact;
    int referenceCount;
    int body1, body2;
    Vec3 location, normal;
    Real radius, depth;
};

/**
 * This is the internal implementation class for TriangleMeshContact.
 */

class SimTK_SIMBODY_EXPORT TriangleMeshContactImpl : public ContactImpl {
public:
    TriangleMeshContactImpl(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth, const std::set<int>& faces1, const std::set<int>& faces2);
private:
    friend class TriangleMeshContact;
    const std::set<int> faces1;
    const std::set<int> faces2;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_IMPL_H_
