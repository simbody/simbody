#ifndef SimTK_SIMBODY_VISUALIZER_H_
#define SimTK_SIMBODY_VISUALIZER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "simbody/internal/common.h"

#include <utility> // for std::pair

namespace SimTK {

class MultibodySystem;
class VisualizationEventListener;
class DecorationGenerator;

class SimTK_SIMBODY_EXPORT Visualizer {
public:
    /** Construct new Visualizer using default window title (executable name). **/
    Visualizer(MultibodySystem& system);
    /** Construct new Visualizer with a given window title. **/
    Visualizer(MultibodySystem& system, const String& title);
    ~Visualizer();
    void report(const State& state) const;
    void addEventListener(VisualizationEventListener* listener);
    const Array_<VisualizationEventListener*>& getEventListeners() const;
    void addMenu(const String& title, const Array_<std::pair<String, int> >& items);
    /**
     * Add an always-present, body-fixed piece of geometry like the one passed in, but attached to the
     * indicated body. The supplied transform is applied on top of whatever transform is already contained
     * in the supplied geometry, and any body Id stored with the geometry is ignored.
     */
    void addDecoration(MobilizedBodyIndex, const Transform& X_BD, const DecorativeGeometry&);

    /**
     * Add an always-present rubber band line, modeled after the DecorativeLine supplied here.
     * The end points of the supplied line are ignored, however: at run time the spatial locations
     * of the two supplied stations are calculated and used as end points.
     */
    void addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1,
                           MobilizedBodyIndex b2, const Vec3& station2,
                           const DecorativeLine& line);
    /**
     * Add a DecorationGenerator that will be invoked to add dynamically generated geometry
     * to the scene.  The Visualizer assumes ownership of the object passed to this method,
     * and will delete it when the Visualizer is deleted.
     */
    void addDecorationGenerator(DecorationGenerator* generator);
    /**
     * Set the transform defining the position and orientation of the camera.
     */
    void setCameraTransform(const Transform& transform);
    /**
     * Move the camera forward or backward so that all geometry in the scene is visible.
     */
    void zoomCameraToShowAllGeometry();
    /**
     * Rotate the camera so that it looks at a specified point.
     *
     * @param point       the point to look at
     * @param upDirection a direction which should point upward as seen by the camera
     */
    void pointCameraAt(const Vec3& point, const Vec3& upDirection);
    /**
     * Set the camera's vertical field of view, measured in radians.
     */
    void setCameraFieldOfView(Real fov);
    /**
     * Set the distance from the camera to the near and far clipping planes.
     */
    void setCameraClippingPlanes(Real nearPlane, Real farPlane);
    /**
     * Set the position and orientation of the ground plane.
     *
     * @param axis     the axis the ground plane is perpendicular to
     * @param height   the position of the ground plane along the specified axis
     */
    void setGroundPosition(const CoordinateAxis& axis, Real height);
    class VisualizerRep;

    // OBSOLETE NAME: will be removed in a later release.
    void zoomCameraToIncludeAllGeometry() {zoomCameraToShowAllGeometry();}
private:
    class RubberBandLine;
    VisualizerRep* rep;
    const VisualizerRep& getRep() const {assert(rep); return *rep;}
    VisualizerRep&       updRep() const {assert(rep); return *rep;}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_VISUALIZER_H_
