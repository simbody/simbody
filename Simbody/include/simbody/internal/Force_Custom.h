#ifndef SimTK_SIMBODY_FORCE_CUSTOM_H_
#define SimTK_SIMBODY_FORCE_CUSTOM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
 * Contributors: Nabeel Allana, Chris Dembia, Thomas Lau                      *
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

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force
subclass Force::Custom, and its related class Force::Custom::Implemnetation,
and is logically part of Force.h. The file assumes that
Force.h will have included all necessary declarations. **/

namespace SimTK {

/**
 * This class is used to define new force elements. To use it, you will
 * create a class that extends Force::Custom::Implementation. Optionally, you
 * may also create a "handle" class extending Force::Custom that hides your
 * implementation and provides a nicer API for users of your new force element.
 * Here we'll use as an example a force "MySpring" that we'll presume you will
 * declare in a header file "MySpring.h". The MySpring constructor will take
 * a reference to a MobilizedBody \a mobod and a spring constant \a k, and
 * connect the origin OB of body \a mobod to the Ground origin O by a spring of
 * stiffness \a k. MySpring will apply a force of magnitude \a kx to OB,
 * directed towards O, where x=|OB-O| is the current distance from OB to O.
 * First we'll look at how MySpring would be used in a main program (whether by
 * you or some future user of your force element), then how you would
 * implement it.
 *
 * There are two possibilities:
 *  - you supply only an %Implementation object MySpringImpl (less work
 *    for you), or
 *  - you supply both the %Implementation object and a handle MySpring
 *    (nicer for the user).
 *
 * If you are planning only to use this force yourself, and perhaps just once
 * for a single application, the %Implementation-only approach is probably
 * adequate. However, if you plan to have others use your new force object, it
 * is well worth the (minimal) extra effort to provide a handle class also to
 * make your new force behave identically to the built-in forces that come
 * with Simbody. In either case the %Implementation object is the same so you
 * can add a handle later if you want. To reiterate: the handle class is
 * completely optional; you \e must write an %Implementation class
 * but a handle class is an aesthetic addition whose main purpose is to
 * make a cleaner API for users of your force element.
 *
 * In the case where you write only the %Implementation class, a user will
 * create an instance of that class and pass it to the generic Force::Custom
 * constructor (this would be in main.cpp or some other piece of application
 * code):
 * @code
 *  // Using your custom force in a simulation application, no handle.
 *  #include "Simbody.h"
 *  #include "MySpring.h"
 *  using namespace SimTK;
 *  GeneralForceSubsystem forces;
 *  MobilizedBody mobod1 = MobilizedBody::Pin(...); // or whatever
 *  // ...
 *  Force::Custom spring1(forces, new MySpringImpl(mobod1,100.));
 * @endcode
 * If you also write a handle class, then the API seen by the user is the
 * same as for built-in force objects:
 * @code
 *  // Using your custom force in a simulation application, with handle.
 *  // Same as above except last line is now:
 *  MySpring spring1(forces, mobod1, 100.);
 * @endcode
 *
 * <h2>Writing the %Implementation class</h2>
 *
 * You would put the following code in MySpring.h, although you can hide it
 * elsewhere (such as in a separate .cpp file) if you provide a handle class:
 * @code
 *  // Sample implementation of a custom force Implementation class.
 *  class MySpringImpl: public Force::Custom::Implementation {
 *  public:
 *      MySpringImpl(MobilizedBody mobod, Real k)
 *      :   m_mobod(mobod), m_k(k) {}
 *
 *      virtual void calcForce(const State&         state,
 *                             Vector_<SpatialVec>& bodyForcesInG,
 *                             Vector_<Vec3>&       particleForcesInG,
 *                             Vector&              mobilityForces) const
 *      {
 *          Vec3 bodyPointInB(0,0,0); // body origin (in the body frame)
 *          Vec3 pos = m_mobod.getBodyOriginLocation(state); // in Ground frame
 *          // apply force towards origin, proportional to distance x
 *          m_mobod.applyForceToBodyPoint(state, bodyPointInB, -m_k*pos,
 *                                        bodyForcesInG);
 *      }
 *
 *      virtual Real calcPotentialEnergy(const State& state) const {
 *          Vec3 pos = m_mobod.getBodyOriginLocation(state); // in Ground
 *          Real x   = pos.norm(); // distance from Ground origin
 *          return m_k*x*x/2; // potential energy in the spring
 *      }
 *  private:
 *      MobilizedBody m_mobod;    // the body to which to apply the force
 *      Real          m_k;        // spring stiffness
 *  };
 * @endcode
 *
 * To write the code exactly as above, the compiler has to be told to look in
 * the %SimTK namespace for some of the symbols. There are a variety of ways
 * to do that; see the discussion below for details.
 *
 * <h2>Writing the handle class</h2>
 *
 * Here is how to implement the handle class, assuming you've already
 * written the %Implementation class MySpringImpl. You would put this
 * code (at least the declarations) in MySpring.h, following the declaration
 * of the MySpringImpl class:
 * @code
 *  // Sample implementation of a handle class. Note: a handle class
 *  // may *not* have any data members or virtual methods.
 *  class MySpring : public Force::Custom {
 *  public:
 *      MySpring(GeneralForceSubsystem& forces, // note the "&"
 *               MobilizedBody          mobod,
 *               Real                   k)
 *      :   Force::Custom(forces, new MySpringImpl(mobod,k)) {}
 *  };
 * @endcode
 * As you can see, the handle class is very simple and just hides the creation
 * of the implementation object. Since this removes any reference to the
 * implementation object from the user's program, it also means you can hide
 * the implementation details completely (perhaps in a separately compiled
 * library), which has many advantages. You can add additional methods to the
 * handle class to provide a clean API to users of your custom force; these
 * methods will forward to the implementation object as necessary but will not
 * expose any aspects of the implementation that are not needed by the user.
 *
 * @warning A handle class <em>must not</em> have any data members or virtual
 * methods and will not work properly if it does. Instead, store any data you
 * need in the %Implementation object.
 *
 * <h2>Getting access to symbols in the %SimTK namespace</h2>
 *
 * The examples above glossed over the naming of symbols in the SimTK namespace.
 * To write the code exactly as above you would need to precede it with
 * \c using statements to tell the compiler to look there to resolve symbols,
 * for example:
 * @code
 *  using namespace SimTK; // search the entire namespace for symbols, or
 *  using SimTK::Real; using SimTK::Vector_; // define just particular symbols
 * @endcode
 * Either of those approaches will work, but note that this will have the same
 * effect on user programs that include MySpring.h as it does within the header
 * file. That may be unwanted behavior for some users who might prefer not to have
 * the namespace searched automatically, perhaps to avoid conflicts with their own
 * symbols that have the same names. If you want to avoid introducing unwanted
 * symbols into users' compilation units, then in the header file you should
 * refer to each symbol explicitly by its full name; that is, prepend the
 * namespace each time the symbol is used, for example:
 * @code
 *  class MySpringImpl: public SimTK::Force::Custom::Implementation {
 *  void calcForce(const SimTK::State&                state,
 *                 SimTK::Vector_<SimTK::SpatialVec>& bodyForcesInG,
 *                 SimTK::Vector_<SimTK::Vec3>&       particleForcesInG,
 *                 SimTK::Vector&                     mobilityForces) const ;
 *  };
 * @endcode
 * That is less convenient for you (and uglier) but avoids the possibility of
 * unwanted side effects in user code so should be considered if you expect
 * wide distribution of your new force element.
 *
 * Thanks to Nur Adila Faruk Senan (a.k.a. Adila Papaya) for help in clarifying
 * this documentation.
 */
class SimTK_SIMBODY_EXPORT Force::Custom : public Force {
public:
    class Implementation;
    /**
     * Create a Custom force.
     *
     * @param forces         the subsystem to which this force should be added
     * @param implementation the object which implements the custom force.  The Force::Custom takes over
     *                       ownership of the implementation object, and deletes it when the Force itself
     *                       is deleted.
     */
    Custom(GeneralForceSubsystem& forces, Implementation* implementation);

    /** Default constructor creates an empty handle. **/
    Custom() {}

    /** @cond **/ // don't let Doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, Force);
    /** @endcond **/
protected:
    const Implementation& getImplementation() const;
    Implementation& updImplementation();
};

/**
 * Every custom force requires implementation of a class that is derived from
 * this abstract class.\ See Force::Custom for details.
 */
class SimTK_SIMBODY_EXPORT Force::Custom::Implementation {
public:
    virtual ~Implementation() { }
    /**
     * Calculate the force for a given state.
     *
     * @param state          the State for which to calculate the force
     * @param bodyForces     spatial forces on MobilizedBodies are accumulated in this.  To apply a force to a body,
     *                       add it to the appropriate element of this vector.
     * @param particleForces forces on particles are accumulated in this.  Since particles are not yet implemented,
     *                       this is ignored.
     * @param mobilityForces forces on individual mobilities (elements of the state's u vector) are accumulated in this.
     *                       To apply a force to a mobility, add it to the appropriate element of this vector.
     */
    virtual void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const = 0;
    /**
     * Calculate this force's contribution to the potential energy of the System.
     *
     * @param state          the State for which to calculate the potential energy
     */
    virtual Real calcPotentialEnergy(const State& state) const = 0;
    /**
     * Get whether this force depends only on the position variables (q), not on the velocies (u) or auxiliary variables (z).
     * The default implementation returns false.  If the force depends only on positions, you should override this to return
     * true.  This allows force calculations to be optimized in some cases.
     */
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
    /**
     * Called during the
     */
    virtual bool isParallelByDefault() const {
        return false;
    }
    /** The following methods may optionally be overridden to do specialized
    realization for a Force. **/
    //@{
    virtual void realizeTopology(State& state) const {}
    virtual void realizeModel(State& state) const {}
    virtual void realizeInstance(const State& state) const {}
    virtual void realizeTime(const State& state) const {}
    virtual void realizePosition(const State& state) const {}
    virtual void realizeVelocity(const State& state) const {}
    virtual void realizeDynamics(const State& state) const {}
    virtual void realizeAcceleration(const State& state) const {}
    virtual void realizeReport(const State& state) const {}
    //@}

    /** Override this if you want to generate some geometry for the visualizer
    to display. Be sure to \e append the geometry to the list. **/
    virtual void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage,
        Array_<DecorativeGeometry>& geometry) const {}

};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_CUSTOM_H_
