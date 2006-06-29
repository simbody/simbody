#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"

#include <vector>

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

class MatterSubsystem;
class ForceSubsystem;


/**
 * The job of the MultibodySystem class is to coordinate the activities of various
 * subsystems which can be part of a multibody system. We insist on having exactly one
 * MatterSubsystem, and we would like also to have:
 *    - a ForceSubsystem
 *    - an AnalyticGeometrySubsystem
 *    - a MassPropertiesSubsystem
 *    - a VisualizationSubsystem
 * There will also be a generic System-level "subsystem" for global variables.
 */
class SimTK_SIMBODY_API MultibodySystem : public System {
public:
    MultibodySystem();
    MultibodySystem(MatterSubsystem& m, ForceSubsystem& f);

    // We inherit realize() from System, and add constraint projection here.
    // We are given a state whose continuous state variables y may violate
    // a set of constraints at position (q) and velocity (u) levels. In addition
    // we may be given a set of absolute error estimates y_err for y. This solver
    // performs two operations:
    //   (1) perform a least squares projection of y onto the constraint manifold,
    //       using the error test norm to define the least-squares direction
    //   (2) perform the same projection on y_err, returning a revised y_err which
    //       has a smaller norm
    // This routine returns true if any change was made to s or y_err, otherwise false.
    // 

    bool project(State& s, Vector& y_err, 
                 const Real& tol,               // must achieve this tolerance or better
                 const Real& dontProjectFac,    // skip projection if tol <= fac*tol
                 const Real& targetTol          // when projecting, try for this (<= tol)
                 ) const;


    // Steals ownership of the source.
    MatterSubsystem& addMatterSubsystem(MatterSubsystem&);
    ForceSubsystem&  addForceSubsystem(ForceSubsystem&);

    int getNMatterSubsystems() const;
    int getNForceSubsystems()  const;
    const MatterSubsystem& getMatterSubsystem(int i) const;
    const ForceSubsystem&  getForceSubsystem(int i)  const;
    MatterSubsystem& updMatterSubsystem(int i);
    ForceSubsystem&  updForceSubsystem(int i);

    // Global state variables dealing with interaction between forces & matter

    // Responses available when the global subsystem is advanced to Dynamics stage.
    const Vector_<SpatialVec>& getRigidBodyForces(const State&,int matterSubsysNum) const;
    const Vector_<Vec3>&       getParticleForces(const State&,int matterSubsysNum) const;
    const Vector&              getMobilityForces(const State&,int matterSubsysNum) const;
    const Real&                getPotentialEnergy(const State&) const;
    const Real&                getKineticEnergy(const State&) const;

    // TODO: camera facing, screen fixed, calculated geometry (e.g. line between stations
    // on two different bodies, marker at system COM)
    void addAnalyticGeometry  (int body, const Transform& X_BG, const AnalyticGeometry&);
    void addDecorativeGeometry(int body, const Transform& X_BG, const DecorativeGeometry&);
    const Array<AnalyticGeometry>&   getBodyAnalyticGeometry(int body);
    const Array<DecorativeGeometry>& getBodyDecorativeGeometry(int body);

    SimTK_PIMPL_DOWNCAST(MultibodySystem, System);
    class MultibodySystemRep& updRep();
    const MultibodySystemRep& getRep() const;
};

class SimTK_SIMBODY_API MultibodyDynamicsStudy : public Study {
public:
    MultibodyDynamicsStudy() { }
    MultibodyDynamicsStudy(const MultibodyDynamicsStudy&);
    MultibodyDynamicsStudy& operator=(const MultibodyDynamicsStudy&);
    ~MultibodyDynamicsStudy();

    MultibodyDynamicsStudy(const MultibodySystem&);

    const MultibodySystem& getMultibodySystem() const;

    void advanceTimeBy(const Real& h); //TODO

    SimTK_PIMPL_DOWNCAST(MultibodyDynamicsStudy, Study);
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
