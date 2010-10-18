#ifndef SimTK_SIMBODY_VISUALIZATION_GEOMETRY_H_
#define SimTK_SIMBODY_VISUALIZATION_GEOMETRY_H_

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

/** @file
 * This is the implementation of DecorativeGeometry used by VisualzationRporter.
 */

#include "SimTKcommon.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Visualizer.h"

namespace SimTK {

class VisualizationGeometry : public DecorativeGeometryImplementation {
public:
    VisualizationGeometry(const Visualizer& visualizer, const SimbodyMatterSubsystem& matter, const State& state);
    ~VisualizationGeometry() {
    }
    void implementLineGeometry(const DecorativeLine& geom);
    void implementBrickGeometry(const DecorativeBrick& geom);
    void implementCylinderGeometry(const DecorativeCylinder& geom);
    void implementCircleGeometry(const DecorativeCircle& geom);
    void implementSphereGeometry(const DecorativeSphere& geom);
    void implementEllipsoidGeometry(const DecorativeEllipsoid& geom);
    void implementFrameGeometry(const DecorativeFrame& geom);
    void implementTextGeometry(const DecorativeText& geom);
    void implementMeshGeometry(const DecorativeMesh& geom);
private:
    Vec4 getColor(const DecorativeGeometry& geom) const;
    int getRepresentation(const DecorativeGeometry& geom) const;
    const Visualizer& visualizer;
    const SimbodyMatterSubsystem& matter;
    const State& state;
};

}

#endif // SimTK_SIMBODY_VISUALIZATION_GEOMETRY_H_
