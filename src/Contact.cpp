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


#include "simbody/internal/ContactImpl.h"
#include <set>

using namespace SimTK;
using std::set;

ContactImpl::ContactImpl(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth) :
        body1(body1), body2(body2), location(location), normal(normal), radius(radius), depth(depth), referenceCount(0) {
}

Contact::Contact(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth) :
        impl(new ContactImpl(body1, body2, location, normal, radius, depth)) {
    impl->referenceCount++;
}

Contact::Contact(ContactImpl* impl) : impl(impl) {
    impl->referenceCount++;
}

Contact::Contact(const Contact& copy) : impl(copy.impl) {
    impl->referenceCount++;
}

Contact& Contact::operator=(const Contact& copy) {
    if (impl) {
        impl->referenceCount--;
        if (impl->referenceCount == 0)
            delete impl;
    }
    impl = copy.impl;
    assert(impl);
    impl->referenceCount++;
    return *this;
}

Contact::~Contact() {
    if (impl) {
        impl->referenceCount--;
        if (impl->referenceCount == 0)
            delete impl;
    }
}

int Contact::getFirstBody() const {
    return impl->body1;
}

int Contact::getSecondBody() const {
    return impl->body2;
}

Vec3 Contact::getLocation() const {
    return impl->location;
}

Vec3 Contact::getNormal() const {
    return impl->normal;
}

Real Contact::getRadius() const {
    return impl->radius;
}

Real Contact::getDepth() const {
    return impl->depth;
}

TriangleMeshContactImpl::TriangleMeshContactImpl(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth,
        const set<int>& faces1, const set<int>& faces2) : ContactImpl(body1, body2, location, normal, radius, depth),
        faces1(faces1), faces2(faces2) {
}

TriangleMeshContact::TriangleMeshContact(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth, const std::set<int>& faces1, const std::set<int>& faces2) :
        Contact(new TriangleMeshContactImpl(body1, body2, location, normal, radius, depth, faces1, faces2)) {
}

const set<int>& TriangleMeshContact::getFirstBodyFaces() const {
    return dynamic_cast<const TriangleMeshContactImpl&>(getImpl()).faces1;
}

const set<int>& TriangleMeshContact::getSecondBodyFaces() const {
    return dynamic_cast<const TriangleMeshContactImpl&>(getImpl()).faces2;
}

bool TriangleMeshContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const TriangleMeshContactImpl*>(&contact.getImpl()) != 0);
}
