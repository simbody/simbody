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

#include "simbody/internal/VisualizationGeometry.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

using namespace SimTK;

static const Vec3 DefaultBodyColor = Gray;

VisualizationGeometry::VisualizationGeometry(const VisualizationProtocol& protocol, const SimbodyMatterSubsystem& matter, const State& state) :
        protocol(protocol), matter(matter), state(state) {
}

void VisualizationGeometry::implementLineGeometry(const SimTK::DecorativeLine& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawLine(transform*geom.getPoint1(), transform*geom.getPoint2(), getColor(geom), geom.getLineThickness() == -1 ? 1 : geom.getLineThickness());
}

void VisualizationGeometry::implementBrickGeometry(const SimTK::DecorativeBrick& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawBox(transform*geom.getTransform(), getScale(geom)*geom.getHalfLengths(), getColor(geom), getRepresentation(geom));
}

void VisualizationGeometry::implementCylinderGeometry(const SimTK::DecorativeCylinder& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawCylinder(transform*geom.getTransform(), getScale(geom)*Vec3(geom.getRadius(), geom.getHalfHeight(), geom.getRadius()), getColor(geom), getRepresentation(geom));
}

void VisualizationGeometry::implementCircleGeometry(const SimTK::DecorativeCircle& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawCircle(transform*geom.getTransform(), getScale(geom)*Vec3(geom.getRadius(), geom.getRadius(), 1), getColor(geom), getRepresentation(geom));
}

void VisualizationGeometry::implementSphereGeometry(const SimTK::DecorativeSphere& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawEllipsoid(transform*geom.getTransform(), getScale(geom)*Vec3(geom.getRadius()), getColor(geom), getRepresentation(geom));
}

void VisualizationGeometry::implementEllipsoidGeometry(const SimTK::DecorativeEllipsoid& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawEllipsoid(transform*geom.getTransform(), getScale(geom)*geom.getRadii(), getColor(geom), getRepresentation(geom));
}

void VisualizationGeometry::implementFrameGeometry(const SimTK::DecorativeFrame& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawFrame(transform*geom.getTransform(), getScale(geom)*geom.getAxisLength(), getColor(geom));
}

void VisualizationGeometry::implementTextGeometry(const SimTK::DecorativeText& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawText(transform*geom.getTransform().T(), getScale(geom), getColor(geom), geom.getText());
}

void VisualizationGeometry::implementMeshGeometry(const SimTK::DecorativeMesh& geom) {
    const Transform& transform  = matter.getMobilizedBody(MobilizedBodyIndex(geom.getBodyId())).getBodyTransform(state);
    protocol.drawPolygonalMesh(geom.getMesh(), transform*geom.getTransform().T(), getScale(geom), getColor(geom), getRepresentation(geom));
}

Vec4 VisualizationGeometry::getColor(const DecorativeGeometry& geom) {
    Vec4 result;
    result.updSubVec<3>(0) = (geom.getColor()[0] == -1 ? DefaultBodyColor : geom.getColor());
    result[3] = (geom.getOpacity() < 0 ? 1 : geom.getOpacity());
    return result;
}

int VisualizationGeometry::getRepresentation(const DecorativeGeometry& geom) const {
    if (geom.getRepresentation() == DecorativeGeometry::DrawDefault)
        return DecorativeGeometry::DrawSurface;
    return geom.getRepresentation();
}

Real VisualizationGeometry::getScale(const DecorativeGeometry& geom) const {
    if (geom.getScale() == -1)
        return 1;
    return geom.getScale();
}
