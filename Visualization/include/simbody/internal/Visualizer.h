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

#include <map>
#include <string>
#include <vector>

namespace SimTK {

class VisualizationEventListener;

class Visualizer {
public:
    Visualizer();
    void beginScene() const;
    void finishScene() const;
    void drawBox(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const;
    void drawEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const;
    void drawCylinder(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const;
    void drawCircle(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const;
    void drawPolygonalMesh(const PolygonalMesh& mesh, const Transform& transform, Real scale, const Vec4& color, int representation) const;
    void drawLine(const Vec3& end1, const Vec3& end2, const Vec4& color, Real thickness) const;
    void drawText(const Vec3& position, Real scale, const Vec4& color, const std::string& string) const;
    void drawFrame(const Transform& transform, Real axisLength, const Vec4& color) const;
    void addEventListener(VisualizationEventListener* listener);
    const std::vector<VisualizationEventListener*>& getEventListeners() const;
private:
    void drawMesh(const Transform& transform, const Vec3& scale, const Vec4& color, short representation, short meshIndex) const;
    int outPipe;
    mutable std::map<const void*, int> meshes;
    std::vector<VisualizationEventListener*> listeners;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_VISUALIZER_H_
