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
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "VisualizerGeometry.h"
#include "VisualizerProtocol.h"

using namespace SimTK;

static const Vec3 DefaultBodyColor = Gray;
static const Vec3 DefaultPointColor = Magenta;

VisualizerGeometry::VisualizerGeometry
   (VisualizerProtocol& protocol, const SimbodyMatterSubsystem& matter, 
    const State& state) 
:   protocol(protocol), matter(matter), state(state) {}

// The DecorativeGeometry's frame D is given in the body frame B, via transform
// X_BD. We want to know X_GD, the pose of the geometry in Ground, which we get 
// via X_GD=X_GB*X_BD.
Transform VisualizerGeometry::calcX_GD(const DecorativeGeometry& geom) const {
    const MobilizedBody& mobod = 
        matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId()));
    const Transform& X_GB  = mobod.getBodyTransform(state);
    const Transform& X_BD  = geom.getTransform();
    return X_GB*X_BD;
}


// We're going to draw three short lines aligned with the body axes
// that intersect at the point.
void VisualizerGeometry::
implementPointGeometry(const SimTK::DecorativePoint& geom) {
    const MobilizedBody& mobod = 
        matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId()));
    const Transform& X_GB  = mobod.getBodyTransform(state);
    const Transform& X_BD  = geom.getTransform();
    const Transform X_GD = X_GB*X_BD;
    const Vec3 p_GP = X_GD*geom.getPoint();
    const Real thickness = 
        geom.getLineThickness() == -1 ? Real(1) : geom.getLineThickness();

    const Real DefaultLength = 0.05; // 1/20 of a unit length
    const Vec3 lengths = DefaultLength * getScaleFactors(geom);
    const Vec4 color = getColor(geom, DefaultPointColor);

    protocol.drawLine(p_GP - lengths[0]*X_GB.x(), 
                      p_GP + lengths[0]*X_GB.x(), color, thickness);
    protocol.drawLine(p_GP - lengths[1]*X_GB.y(), 
                      p_GP + lengths[1]*X_GB.y(), color, thickness);
    protocol.drawLine(p_GP - lengths[2]*X_GB.z(), 
                      p_GP + lengths[2]*X_GB.z(), color, thickness);
}

void VisualizerGeometry::implementLineGeometry(const SimTK::DecorativeLine& geom) {
    const Transform X_GD = calcX_GD(geom);
    protocol.drawLine(X_GD*geom.getPoint1(), X_GD*geom.getPoint2(), getColor(geom), geom.getLineThickness() == -1 ? 1 : geom.getLineThickness());
}

void VisualizerGeometry::implementBrickGeometry(const SimTK::DecorativeBrick& geom) {
    const Transform X_GD = calcX_GD(geom);
    const Vec3 hlen = getScaleFactors(geom).elementwiseMultiply(geom.getHalfLengths());
    protocol.drawBox(X_GD, hlen, getColor(geom), getRepresentation(geom));
}

void VisualizerGeometry::implementCylinderGeometry(const SimTK::DecorativeCylinder& geom) {
    const Transform X_GD = calcX_GD(geom);
    const Vec3 scale = getScaleFactors(geom);
    protocol.drawCylinder(X_GD, Vec3(scale[0]*geom.getRadius(), 
                                     scale[1]*geom.getHalfHeight(), 
                                     scale[2]*geom.getRadius()), 
                          getColor(geom), getRepresentation(geom), getResolution(geom));
}

void VisualizerGeometry::implementCircleGeometry(const SimTK::DecorativeCircle& geom) {
    const Transform X_GD = calcX_GD(geom);
    const Vec3 scale = getScaleFactors(geom); // z ignored
    protocol.drawCircle(X_GD, Vec3(scale[0]*geom.getRadius(), 
                                   scale[1]*geom.getRadius(), 1), 
                        getColor(geom), getRepresentation(geom), getResolution(geom));
}

void VisualizerGeometry::implementSphereGeometry(const SimTK::DecorativeSphere& geom) {
    const Transform X_GD = calcX_GD(geom);
    protocol.drawEllipsoid(X_GD, geom.getRadius()*getScaleFactors(geom), 
                           getColor(geom), getRepresentation(geom), getResolution(geom));
}

void VisualizerGeometry::implementEllipsoidGeometry(const SimTK::DecorativeEllipsoid& geom) {
    const Transform X_GD = calcX_GD(geom);
    const Vec3 radii = getScaleFactors(geom).elementwiseMultiply(geom.getRadii());
    protocol.drawEllipsoid(X_GD, radii, getColor(geom), getRepresentation(geom),
                           getResolution(geom));
}

void VisualizerGeometry::implementFrameGeometry(const SimTK::DecorativeFrame& geom) {
    const Transform X_GD = calcX_GD(geom);
    protocol.drawCoords(X_GD, geom.getAxisLength()*getScaleFactors(geom), getColor(geom));
}

void VisualizerGeometry::implementTextGeometry(const SimTK::DecorativeText& geom) {
    const Transform X_GD = calcX_GD(geom);
    // The default is to face the camera.
    bool faceCamera = geom.getFaceCamera()<0 ? true : (geom.getFaceCamera()!=0);
    protocol.drawText(X_GD.p(), getScaleFactors(geom), getColor(geom), 
                      geom.getText(), faceCamera);
}

void VisualizerGeometry::implementMeshGeometry(const SimTK::DecorativeMesh& geom) {
    const Transform X_GD = calcX_GD(geom);
    protocol.drawPolygonalMesh(geom.getMesh(), X_GD, getScaleFactors(geom), 
                               getColor(geom), getRepresentation(geom));
}

Vec4 VisualizerGeometry::getColor(const DecorativeGeometry& geom,
                                  const Vec3& defaultColor) {
    Vec4 result;
    if (geom.getColor()[0] >= 0) 
        result.updSubVec<3>(0) = geom.getColor();
    else {
        const Vec3 def = defaultColor[0] >= 0 ? defaultColor : DefaultBodyColor;
        result.updSubVec<3>(0) = def;
    }
    result[3] = (geom.getOpacity() < 0 ? 1 : geom.getOpacity());
    return result;
}

int VisualizerGeometry::getRepresentation(const DecorativeGeometry& geom) const {
    if (geom.getRepresentation() == DecorativeGeometry::DrawDefault)
        return DecorativeGeometry::DrawSurface;
    return geom.getRepresentation();
}

unsigned short VisualizerGeometry::
getResolution(const DecorativeGeometry& geom) const {
    if (geom.getResolution() <= 0)
        return 2;
    return std::max((unsigned short) 1, (unsigned short) (geom.getResolution()*2));
}

Vec3 VisualizerGeometry::getScaleFactors(const DecorativeGeometry& geom) const {
    const Vec3& scale = geom.getScaleFactors();
    Vec3 actual;
    for (int i=0; i<3; ++i)
        actual[i] = scale[i] <= 0 ? 1 : scale[i];
    return actual;
}
