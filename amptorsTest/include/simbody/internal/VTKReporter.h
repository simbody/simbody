#ifndef SimTK_SIMBODY_VTK_REPORTER_H_
#define SimTK_SIMBODY_VTK_REPORTER_H_

/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/** @file
 * This is the user-visible handle class for the VTK Reporter which
 * provides a 3d animation window for viewing MultibodySystem
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

class SimTK_SIMBODY_EXPORT VTKReporter  {
public:
    VTKReporter() : rep(0) { }

    // Set the scale to 0 to disable automatically-generated geometry. Otherwise set
    // it to a typical length scale for a body.
    explicit VTKReporter(const MultibodySystem& m, Real defaultScaleForAutoGeometry=1.);

    VTKReporter(const VTKReporter&);
    ~VTKReporter();
    VTKReporter& operator=(const VTKReporter&);

    /// This method calculates a new visualization frame using the MultibodySystem
    /// and the supplied state, and updates the screen.
    void report(const State& s);

    /// Change the camera location. This is a point in the Ground frame.
    void setCameraLocation(const Vec3&);
    /// Change the point towards which the camera is pointing. This point is
    /// in the Ground frame.
    void setCameraFocalPoint(const Vec3&);
    /// Change the camera's "up" direction. This is a vector in the ground frame;
    /// only its direction is meaningful.
    void setCameraUpDirection(const Vec3&);
    /// Set the location of the clipping planes. These planes are perpendicular to the
    /// view plane normal, which is a vector in the direction from the camera location
    /// to the camera focal point. Their locations are given by distances along the
    /// the view plane normal, measured from the camera position.
    void setCameraClippingRange(Real nearPlane, Real farPlane);
    /// Zoom the camera so that all geometry is visible.
    void zoomCameraToIncludeAllGeometry();
    /// Zoom by an amount relative to the current zoom.
    void zoomCamera(Real);

    /// Add an always-present, body-fixed piece of geometry like the one passed in, but attached to the
    /// indicated body. The supplied transform is applied on top of whatever transform is already contained
    /// in the supplied geometry, and any body Id stored with the geometry is ignored.
    /// The 3d representation of the geometry here can be precalculated; only the orientation
    /// of the body frame needs to be applied at run time.
    void addDecoration(BodyId, const Transform& X_BD, const DecorativeGeometry&);

    /// Add an always-present rubber band line, modeled after the DecorativeLine supplied here.
    /// The end points of the supplied line are ignored, however -- at run time we'll calculate
    /// the spatial locations of the two supplied stations and use those as end points. Note
    /// that the 3d representation of this line can't be precalculated because the line length
    /// will vary.
    void addRubberBandLine(BodyId b1, const Vec3& station1, 
                           BodyId b2, const Vec3& station2,
                           const DecorativeLine&);

    /// Add a piece of geometry to the next frame only. The bodies, transforms, etc. are
    /// pulled from the argument, and the geometry is generated immediately and added to
    /// the current frame. Then it is discarded.
    void addEphemeralDecoration(const DecorativeGeometry&);

    // TODO: default geometry generation should be moved to the matter subsystem.
    void disableDefaultGeometry();
    void setDefaultBodyColor(BodyId bodyNum, const Vec3& rgb);
    const Vec3& getDefaultBodyColor(BodyId bodyNum) const;
 
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
