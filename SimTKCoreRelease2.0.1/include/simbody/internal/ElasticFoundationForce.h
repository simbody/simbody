#ifndef SimTK_SIMBODY_ELASTIC_FOUNDATION_FORCE_H_
#define SimTK_SIMBODY_ELASTIC_FOUNDATION_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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


#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/Force.h"

namespace SimTK {

class GeneralContactSubsystem;
class ElasticFoundationForceImpl;

/**
 * This class implements an elastic foundation or "bed of springs" contact model.  It places
 * springs over the surface of a triangle mesh, each of which independently interacts with
 * objects in contact with the mesh.  The interaction includes components for the normal
 * restoring force, dissipation in the material, and surface friction.  This force
 * is only applied to contacts involving a triangle mesh, though the object colliding with
 * the mesh may be of any type.  Contacts which do not involve a triangle mesh are ignored.
 *
 * This class relies on a GeneralContactSubsystem to identify contacts, then applies
 * forces to all contacts in a single contact set.  To use it, do the following:
 *
 * <ol>
 * <li>Add a GeneralContactSubsystem and GeneralForceSubsystem to a MultibodySystem.</li>
 * <li>Create a contact set within the GeneralContactSubsystem, and add ContactGeometry objects to it.</li>
 * <li>Add an ElasticFoundationForce to the GeneralForceSubsystem, and call setBodyParameters() on it
 * once for each TriangleMesh in the contact set.</li>
 * </ol>
 *
 * A spring is placed at the center of each face of each triangle mesh.  When a mesh collides with
 * any other object, the other object is considered to be rigid and the displacement of each spring is
 * calculated independently.  A spring is considered to be displaced if its base location is inside
 * the other object, and the contact point is taken to be the nearest point on that object's surface.
 * When two meshes collide, the springs on each mesh are treated independently: each mesh is assumed
 * to be rigid for purposes of calculating the force exerted by the other mesh's springs.
 *
 * <h1>Normal Force Components</h1>
 *
 * The force exerted by each spring along its displacement direction is given by
 *
 * f = k*a*x*(1+c*v)
 * 
 * where k is the spring stiffness, a is the area of the face the spring belongs to, x is the displacement
 * distance, c is the spring's dissipation coefficient, and v=dx/dt.  If the springs are assumed
 * to represent a uniform layer of elastic material over a rigid substrate, the stiffness is given by
 *
 * k = (1-p)*E/((1+p)(1-2p)*h)
 *
 * where E is the Young's modulus of the elastic layer, p is its Poisson's ratio, and h is its thickness.
 * 
 * <h1>Friction Force</h1>
 *
 * The friction force exerted by each spring is based on a model by Michael Hollars:
 *
 * f = fn*[min(vs/vt,1)*(ud+2(us-ud)/(1+(vs/vt)^2))+uv*vs]
 *
 * where fn is the normal force at the contact point, vs is the slip (tangential)
 * velocity of the two bodies at the contact point, vt is a transition velocity
 * (see below), and us, ud, and uv are the coefficients of static, dynamic, and
 * viscous friction respectively.
 *
 * Because the friction force is a continuous function of the slip velocity, this
 * model cannot represent stiction; as long as a tangential force is applied, the
 * two bodies will move relative to each other.  There will always be a nonzero
 * drift, no matter how small the force is.  The transition velocity vt acts as an
 * upper limit on the drift velocity.  By setting vt to a sufficiently small value,
 * the drift velocity can be made arbitrarily small, at the cost of making the
 * equations of motion very stiff.  The default value of vt is 0.01.
 */
class SimTK_SIMBODY_EXPORT ElasticFoundationForce : public Force {
public:
    /**
     * Create an elastic foundation contact model.
     * 
     * @param forces         the subsystem which will own this ElasticFoundationForce element
     * @param contacts       the subsystem to which this contact model should be applied
     * @param contactSet     the index of the contact set to which this contact model will be applied
     */
    ElasticFoundationForce(GeneralForceSubsystem& forces, GeneralContactSubsystem& contacts, ContactSetIndex contactSet);
    /**
     * Set the material parameters for a body in the contact set, which must be a ContactGeometry::TriangleMesh.
     *
     * @param bodyIndex       the index of the body within the contact set
     * @param stiffness       the stiffness constant (k) for the body
     * @param dissipation     the dissipation coefficient (c) for the body
     * @param staticFriction  the coefficient of static friction (us) for the body
     * @param dynamicFriction the coefficient of dynamic friction (ud) for the body
     * @param viscousFriction the coefficient of viscous friction (uv) for the body
     */
    void setBodyParameters(int bodyIndex, Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction);
    /**
     * Get the transition velocity (vt) of the friction model.
     */
    Real getTransitionVelocity() const;
    /**
     * Set the transition velocity (vt) of the friction model.
     */
    void setTransitionVelocity(Real v);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(ElasticFoundationForce, ElasticFoundationForceImpl, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ELASTIC_FOUNDATION_FORCE_H_
