#ifndef SimTK_SIMBODY_VISUALIZER_GEOMETRY_H_
#define SimTK_SIMBODY_VISUALIZER_GEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Ayman Habib                                                  *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
 * This is the implementation of DecorativeGeometry used by 
 * Visualizer::Reporter.
 */

#include "simbody/internal/common.h"

namespace SimTK {
class SimbodyMatterSubsystem;
class VisualizerProtocol;

class VisualizerGeometry : public DecorativeGeometryImplementation {
public:
    VisualizerGeometry(VisualizerProtocol& protocol, const SimbodyMatterSubsystem& matter, const State& state);
    ~VisualizerGeometry() {
    }
    void implementPointGeometry(const DecorativePoint& geom) override;
    void implementLineGeometry(const DecorativeLine& geom) override;
    void implementBrickGeometry(const DecorativeBrick& geom) override;
    void implementCylinderGeometry(const DecorativeCylinder& geom) override;
    void implementCircleGeometry(const DecorativeCircle& geom) override;
    void implementSphereGeometry(const DecorativeSphere& geom) override;
    void implementEllipsoidGeometry(const DecorativeEllipsoid& geom) override;
    void implementFrameGeometry(const DecorativeFrame& geom) override;
    void implementTextGeometry(const DecorativeText& geom) override;
    void implementMeshGeometry(const DecorativeMesh& geom) override;
    void implementMeshFileGeometry(const DecorativeMeshFile& geom) override {}; // Not handled yet by this Visualizer
    void implementArrowGeometry(const DecorativeArrow& geom) override {}; // Not handled yet by this Visualizer
    void implementTorusGeometry(const DecorativeTorus& geom) override {}; // Not handled yet by this Visualizer
    void implementConeGeometry(const DecorativeCone& geom) override {}; // Not handled yet by this Visualizer
    static Vec4 getColor(const DecorativeGeometry& geom,
                         const Vec3& defaultColor = Vec3(-1));
private:
    int getRepresentation(const DecorativeGeometry& geom) const;
    unsigned short getResolution(const DecorativeGeometry& geom) const;
    Vec3 getScaleFactors(const DecorativeGeometry& geom) const;
    Transform calcX_GD(const DecorativeGeometry& geom) const;
    VisualizerProtocol& protocol;
    const SimbodyMatterSubsystem& matter;
    const State& state;
};
}

#endif // SimTK_SIMBODY_VISUALIZER_GEOMETRY_H_
