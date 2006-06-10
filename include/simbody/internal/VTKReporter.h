#ifndef SimTK_SIMBODY_VTK_REPORTER_H_
#define SimTK_SIMBODY_VTK_REPORTER_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/** @file
 * This is the user-visible handle class for the VTK Reporter which
 * provides a 3d animation window for viewing Mechanical System
 * simulations.
 */

#include "simbody/internal/common.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>


namespace SimTK {

class State;
class MultibodySystem;
class DecorativeGeometry;
class VTKDecoration;

class SimTK_SIMBODY_API VTKReporter {
public:
    VTKReporter() : rep(0) { }
    explicit VTKReporter(const MultibodySystem& m);
    VTKReporter(const VTKReporter&);
    ~VTKReporter();
    VTKReporter& operator=(const VTKReporter&);

    void report(const State& s);

    void addDecoration(int bodyNum, const Transform& X_GD, const DecorativeGeometry&);
    void addRubberBandLine(int b1, const Vec3& station1, int b2, const Vec3& station2,
                           const DecorativeLine&);


    void setDefaultBodyColor(int bodyNum, const Vec3& rgb);
    const Vec3& getDefaultBodyColor(int bodyNum) const;
 
    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit VTKReporter(class VTKReporterRep* r) : rep(r) { }
    bool                  hasRep() const {return rep!=0;}
    const VTKReporterRep& getRep() const {assert(rep); return *rep;}
    VTKReporterRep&       updRep() const {assert(rep); return *rep;}
protected:
    class VTKReporterRep* rep;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_VTK_REPORTER_H_
