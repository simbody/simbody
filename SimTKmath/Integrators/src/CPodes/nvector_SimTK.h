#ifndef SimTK_CPODES_NVECTOR_SimTK_H_
#define SimTK_CPODES_NVECTOR_SimTK_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK CPodes                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/** @file
 * This header defines a Sundials N_Vector implementation which uses
 * SimTK's Vector class internally.
 */

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include "sundials/sundials_nvector.h"

#include <cassert>

using SimTK::Vector;

// This is the object to which the generic N_Vector struct's "content"
// entry points.
class N_VectorContent_SimTK {
public:
    N_VectorContent_SimTK()
      : treatAsConst(false), ownVector(true), data(new Vector()) { 
    }

    N_VectorContent_SimTK(const Vector& v) 
      : treatAsConst(true), ownVector(false),
        data(const_cast<Vector*>(&v)) {
    }

    N_VectorContent_SimTK(Vector& v) 
      : treatAsConst(false), ownVector(false), data(&v) {
    }

    // Copy constructor makes a deep (new) copy.
    N_VectorContent_SimTK(const N_VectorContent_SimTK& nv) 
      : treatAsConst(false), ownVector(true), data(new Vector(*nv.data)) {
    }

    // Assignment fails if target is const, otherwise reallocates as needed.
    N_VectorContent_SimTK& operator=(const N_VectorContent_SimTK& nv) {
        if (&nv != this) {
            assert(!treatAsConst);
            *data = *nv.data;
        }
        return *this;
    }

    ~N_VectorContent_SimTK() {
        if (ownVector) delete data;
        data = 0;
    }

    const Vector& getVector() const {
        assert(data);
        return *data;
    }

    Vector& updVector() {
        assert(data);
        assert(!treatAsConst);
        return *data;
    }

private:
    bool treatAsConst;  // is this logically const?
    bool ownVector;     // if true, destruct Vector along with this
    Vector* data;

};

// This singleton class defines the operations needed by the abstract
// N_Vector in terms of operations performed on SimTK Vectors.
// The "virtual function table" is filled in once and used
// in *all* N_Vector_SimTK objects.
class N_Vector_Ops_SimTK : public _generic_N_Vector_Ops {
public:
    // Return a reference to the singleton object. You can use
    // the address of this reference to check whether a given N_Vector
    // is really an N_Vector_SimTK.
    static const _generic_N_Vector_Ops& getOps() {
        return Ops;
    }

    // Sundials isn't const correct so it is going to want a
    // non-const pointer to the Ops class here, which is const.
    // We'll just have to trust it.
    static _generic_N_Vector_Ops* getSundialsOpsPtr() {
        return const_cast<N_Vector_Ops_SimTK*>(&Ops);
    }

private:
    N_Vector_Ops_SimTK(); // initialize
    SUNDIALS_EXPORT static const N_Vector_Ops_SimTK Ops;
};

// This is the object to which an N_Vector points when we are using
// SimTK's N_Vector implementation.
class N_Vector_SimTK : public _generic_N_Vector {
public:
    // Default constructor creates an owned, empty, resizable vector.
    N_Vector_SimTK() {
        content = (void*)new N_VectorContent_SimTK();
        ops = N_Vector_Ops_SimTK::getSundialsOpsPtr();
    }

    // This constructor produces a read-only *reference* to an existing Vector.
    // Destruction of the N_Vector_SimTK leaves the original Vector as it was.
    N_Vector_SimTK(const Vector& v) {
       content = (void*)new N_VectorContent_SimTK(v);
       ops = N_Vector_Ops_SimTK::getSundialsOpsPtr();
    }

    // This constructor produces a writable *reference* to an existing
    // Vector. Destruction of the N_Vector_SimTK leaves the original
    // Vector in existence, although its value may have changed.
    N_Vector_SimTK(Vector& v) {
       content = (void*)new N_VectorContent_SimTK(v);
       ops = N_Vector_Ops_SimTK::getSundialsOpsPtr();
    }

    // Copy constructor makes a deep copy. The result will be an owned,
    // resizable vector that will evaporate upon destruction of the
    // N_Vector_SimTK object.
    N_Vector_SimTK(const N_Vector_SimTK& nv) {
        assert(nv.ops == &N_Vector_Ops_SimTK::getOps());
        content = (void*)new N_VectorContent_SimTK(nv.getContent());
        ops = N_Vector_Ops_SimTK::getSundialsOpsPtr();
    }

    // Assignment will fail if this isn't writable.
    N_Vector_SimTK& operator=(const N_Vector_SimTK& nv) {
        assert(   ops    == &N_Vector_Ops_SimTK::getOps() 
               && nv.ops == &N_Vector_Ops_SimTK::getOps());
        if (&nv != this)
            updContent() = nv.getContent();
        return *this;
    }

    ~N_Vector_SimTK() {
        assert(ops == &N_Vector_Ops_SimTK::getOps());
        delete reinterpret_cast<N_VectorContent_SimTK*>(content);
        content = 0;
        ops = 0;
    }

    static bool isA(const N_Vector nv) {
        return nv && (nv->ops == &N_Vector_Ops_SimTK::getOps());
    }

    static const N_Vector_SimTK* downcast(const N_Vector nv) {
        assert(isA(nv));
        return reinterpret_cast<const N_Vector_SimTK*>(nv);
    }

    static N_Vector_SimTK* updDowncast(N_Vector nv) {
        assert(isA(nv));
        return reinterpret_cast<N_Vector_SimTK*>(nv);
    }

    static N_Vector_SimTK* nvclone(const N_Vector nv) {
        assert(isA(nv));
        const N_Vector_SimTK& nvs = *reinterpret_cast<const N_Vector_SimTK*>(nv);
        return new N_Vector_SimTK(nvs);
    }

    static const Vector& getVector(const N_Vector nv) {
        assert(isA(nv));
        const N_Vector_SimTK& nvs = *reinterpret_cast<const N_Vector_SimTK*>(nv);
        return nvs.getContent().getVector();
    }

    static Vector& updVector(N_Vector nv) {
        assert(isA(nv));
        N_Vector_SimTK& nvs = *reinterpret_cast<N_Vector_SimTK*>(nv);
        return nvs.updContent().updVector();
    }

private:
    const N_VectorContent_SimTK& getContent() const {
        return *reinterpret_cast<const N_VectorContent_SimTK*>(content);
    }
    N_VectorContent_SimTK& updContent() {
        return *reinterpret_cast<N_VectorContent_SimTK*>(content);
    }

};

#endif // SimTK_NVECTOR_SimTK_H_
