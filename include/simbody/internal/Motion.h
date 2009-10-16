#ifndef SimTK_SIMBODY_MOTION_H_
#define SimTK_SIMBODY_MOTION_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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
 * This defines the Motion class, which is used to specify how the mobilities
 * associated with a particular mobilizer are to be calculated.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {

class SimbodyMatterSubsystem;
class MobilizedBody;
class Motion;
class MotionImpl;

// We only want the template instantiation to occur once. This symbol is defined 
// in the Simbody compilation unit that defines the Motion class but should not 
// be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_MOTION
    extern template class PIMPLHandle<Motion, MotionImpl, true>;
#endif

/**
 * A Motion object belongs to a particular mobilizer and specifies how the 
 * associated motion is to be calculated. There are two independent aspects:
 * (1) at what level is the motion driven (low to high: acceleration, velocity, 
 * or position), and (2) how is the motion at that level specified. Levels lower 
 * than the driven level are also driven; higher levels are free and determined 
 * by integration. This table shows the possibilities:
 * <pre>
 *      Level How driven    Acceleration    Velocity      Position
 *      ----- ----------    ------------    ----------    ----------
 *       Acc  Zero               0           discrete       free
 *        "   Discrete        discrete         free           "
 *        "   Prescribed      a(t,q,u)          "             "
 *        "   Free           from forces        "             "
 *
 *       Vel  Zero               0              0          discrete
 *        "   Discrete           0           discrete        free
 *        "   Prescribed       dv/dt          v(t,q)          "
 *        "   Fast               0           relax(v)         "
 *
 *       Pos  Zero               0              0             0
 *        "   Discrete           0              0          discrete
 *        "   Prescribed      d2p/dt2         dp/dt          p(t)
 *        "   Fast               0              0          relax(p)
 * </pre>
 * There are two duplicates in the above table: specifying acceleration as
 * Zero is the same as specifying velocity as Discrete and specifying
 * velocity as Zero is the same as specifying position as Discrete.
 *
 * For mobilizers with more than one mobility, the associated Motion controls
 * \e all the mobilities and moreover they are all driven at the same level and
 * by the same method. TODO?: provide Motion::Composite for building up a Motion
 * from individual axis specifications.
 *
 * Motion is a PIMPL-style abstract base class, with concrete classes defined
 * for each kind of Motion. There is a set of built-in Motions and a generic 
 * "Custom" Motion (an abstract base class) from which advanced users may derive 
 * their own Motion objects.
 */
    
class SimTK_SIMBODY_EXPORT Motion : public PIMPLHandle<Motion, MotionImpl, true> {
public:
    /// What is the highest level of motion that is driven? Lower levels are
    /// also driven; higher levels are determined by integration.
    enum Level {
        Acceleration = 0, ///< we know udot; integrate to get u, q
        Velocity     = 1, ///< we know u and udot; integrate to get q
        Position     = 2  ///< we know q, u, and udot
    };
    /// Returns a human-readable name corresponding to the given Level; useful
    /// for debugging. If the Level is unrecognized the method will return
    /// some text to that effect rather than crashing.
    static const char* nameOfLevel(Level);

    /// There are several ways to specify the motion at this Level, and the
    /// selected method also determines lower-level motions. Regular is only
    /// permitted when Level==Acceleration, and Fast is not allowed for that
    /// Level.
    enum Method {
        Zero        = 0,
        Discrete    = 1, ///< motion is "slow"; lower levels are zero
        Prescribed  = 2, ///< motion is function of time and state; <level is derivative
        Free        = 3, ///< accel. calculated from forces, pos and vel integrated
        Fast        = 4  ///< motion is "fast"; lower levels are zero
    };
    /// Returns a human-readable name corresponding to the given Method; useful
    /// for debugging. If the Method is unrecognized the method will return
    /// some text to that effect rather than crashing.
    static const char* nameOfMethod(Method);

    Motion() { }
    explicit Motion(MotionImpl* r) : HandleBase(r) { }

	/// Get the MobilizedBody to which this Motion belongs.
	const MobilizedBody& getMobilizedBody() const;

    Level  getLevel(const State&) const;
    Method getLevelMethod(const State&) const;

    // This implements the above table. Given the (level,method) Motion specification,
    // it reports the actual method to be used for each of the three levels.
    void calcAllMethods(const State& s, Method& qMethod, Method& uMethod, Method& udotMethod) const {
        const Level  level       = getLevel(s);
        const Method levelMethod = getLevelMethod(s);
        Method method[3]; // acc, vel, pos
        method[level] = levelMethod;

        switch (level) {
          case Position:
            method[Velocity]=method[Acceleration]= 
                (levelMethod==Prescribed ? Prescribed : Zero);
            break;
          case Velocity:
            method[Acceleration] = (levelMethod==Prescribed ? Prescribed : Zero);
            method[Position]     = (levelMethod==Zero       ? Discrete   : Free);
            break;
          case Acceleration:
            method[Velocity] = (levelMethod==Zero ? Discrete : Free);
            method[Position] = Free;
            break;
          default:
            assert(!"unrecognized level");
        }

        qMethod    = method[Position];
        uMethod    = method[Velocity];
        udotMethod = method[Acceleration];
    }
    
    class Steady;
    class Linear;
    class Polynomial;
    class Composite;
    class Custom;
    
    class SteadyImpl;
    class LinearImpl;
    class PolynomialImpl;
    class CompositeImpl;
    class CustomImpl;
};

/**
 * This non-holonomic Motion object imposes a constant rate on all mobilities.
 */
class SimTK_SIMBODY_EXPORT Motion::Steady : public Motion {
public:
    /**
     * Create a Motion::Steady where all mobilities have the same velocity.
     * 
     * @param[in,out] mobod the MobilizedBody to which this Motion should be added
     * @param[in]     u     the rate to be applied to all mobilities
     */
    Steady(MobilizedBody& mobod, Real u);
    /**
     * Create a Motion::Steady with different velocities for each mobility
     * specified. Any unspecified mobilities will get zero velocity.
     * 
     * @param[in,out] mobod the MobilizedBody to which this Motion should be added
     * @param[in]     u     the rates to be applied to the first N mobilities; the
     *                      rest are set to zero
     */
    template <int N> SimTK_SIMBODY_EXPORT 
    Steady(MobilizedBody& mobod, const Vec<N>& u); // instantiated in library

    Steady& setDefaultRate(Real u);
    Steady& setOneDefaultRate(UIndex, Real u);
    template <int N> SimTK_SIMBODY_EXPORT 
    Steady& setDefaultRates(const Vec<N>& u); // instantiated in library

    Real getDefaultRate(UIndex=UIndex(0)) const;

    void setRate(State&, Real u) const; // all axes set to u
    void setOneRate(State&, UIndex, Real u) const;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Steady, SteadyImpl, Motion);
};


/**
 * This class can be used to define new motions. To use it, create a class that 
 * extends Motion::Custom::Implementation. You can then create an instance of it 
 * and pass it to the Motion::Custom constructor:
 * 
 * <pre>
 * Motion::Custom myMotion(mobod, new MyMotionImplementation());
 * </pre>
 * 
 * Alternatively, you can create a subclass of Motion::Custom which creates the 
 * Implementation itself:
 * 
 * <pre>
 * class MyMotion : public Motion::Custom {
 * public:
 *   MyMotion(MobilizedBody& mobod) 
 *     : Motion::Custom(mobod, new MyMotionImplementation()) {}
 * };
 * </pre>
 * 
 * This allows a user to simply write
 * 
 * <pre>
 * MyMotion(mobod);
 * </pre>
 * 
 * and not worry about implementation classes or creating objects on the heap. If 
 * you do this, your Motion::Custom handle subclass must not have any data members 
 * or virtual methods.  If it does, it will not work correctly. Instead, store all 
 * data in the Implementation subclass.
 */

class SimTK_SIMBODY_EXPORT Motion::Custom : public Motion {
public:
    class Implementation;
    /**
     * Create a Custom motion.
     * 
     * @param mobod          the MobilizedBody to which this Motion should be added
     * @param implementation the object which implements the custom Motion. The 
     *                       Motion::Custom takes over ownership of the 
     *                       implementation object, and deletes it when the Motion 
     *                       itself is deleted.
     */
    Custom(MobilizedBody& mobod, Implementation* implementation);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, Motion);
protected:
    const Implementation& getImplementation() const;
    Implementation& updImplementation();
};

/**
 * This is the abstract base class for Custom Motion implementations.
 */
class SimTK_SIMBODY_EXPORT Motion::Custom::Implementation {
public:
    /// Destructor is virtual.
    virtual ~Implementation() { }

    /// Override this if you want your Motion objects to be copyable.
    virtual Implementation* clone() const {
        SimTK_ERRCHK_ALWAYS(!"unimplemented",
            "Motion::Custom::Implementation::clone()",
            "Concrete Implementation did not supply a clone() method, "
            "but a copy operation was attempted.");
        /*NOTREACHED*/
        return 0;
    }

    /// A Motion prescribes either position, velocity, or acceleration.
    /// When velocity is prescribed, acceleration must also be 
    /// prescribed as the time derivative of the velocity. And, when
    /// position is prescribed, velocity must also be prescribed as
    /// the time derivative of the position (and acceleration as above).
    /// Thus acceleration is \e always prescribed.
    /// Anything not prescribed will be determined by numerical
    /// integration, by relaxation, or by discrete changes driven by
    /// events, depending on whether the associated mobilizer is 
    /// "regular", "fast", or "slow", respectively.
    virtual Motion::Level getLevel(const State&) const = 0;

    virtual Motion::Method getLevelMethod(const State&) const {
        return Motion::Prescribed;
    }

    /// @name Position (Holonomic) prescribed motion virtuals
    ///
    /// These must be defined if the motion method is "Prescribed" and the 
    /// motion level is "Position". In that case q=q(t), qdot, and qdotdot are 
    /// all required. Note that Simbody passes in the number of q's being 
    /// prescribed; make sure you are seeing what you expect.
    //@{
    /// This operator is called during the MatterSubsystem's realize(Time) 
    /// computation. This Motion's own realizeTime() method will already have
    /// been called. The result must depend only on time and earlier-stage
    /// state variables.
    virtual void calcPrescribedPosition(const State& s, int nq, Real* q) const;

    /// Calculate the time derivative of the prescribed positions. The qdots 
    /// calculated here must be the exact time derivatives of the q's returned 
    /// by calcPrescribedPosition(). So the calculation must be limited to 
    /// the same dependencies, plus the current value of this mobilizer's q's
    /// (or the cross-mobilizer transform X_FM because that depends only on 
    /// those q's). Note that we are return qdots, not u's; they are not 
    /// always the same. Simbody knows how to map from qdots to u's when necessary. 
    ///
    /// This operator is called during the MatterSubsystem's realize(Position) 
    /// computation. This Motion's own realizePosition() method will already 
    /// have been called.  
    virtual void calcPrescribedPositionDot(const State& s, int nq, Real* qdot) const;

    /// Calculate the 2nd time derivative of the prescribed positions. The 
    /// qdotdots calculated here must be the exact time derivatives of the qdots 
    /// returned by calcPrescribedPositionDot(). So the calculation must be 
    /// limited to the same dependencies, plus the current value of this 
    /// mobilizer's qdots (or the cross-mobilizer velocity V_FM because that 
    /// depends only on those qdots). Note that we are return qdotdots, not 
    /// udots; they are not always the same. Simbody knows how to map from 
    /// qdotdots to udots when necessary. 
    ///
    /// This operator is called during the MatterSubsystem's realize(Dynamics) 
    /// computation. This Motion's own realizeDynamics() method will already 
    /// have been called.  
    virtual void calcPrescribedPositionDotDot(const State& s, int nq, Real* qdotdot) const;
    //@}

    /// @name Velocity (Nonholonomic) prescribed motion virtuals
    ///
    /// These must be defined if the motion method is "Prescribed" and the 
    /// motion level is "Velocity". In that case u=u(t,q), and udot are 
    /// both required. Note that Simbody passes in the number of u's being 
    /// prescribed; make sure you are seeing what you expect.
    //@{
    /// This operator is called during the MatterSubsystem's realize(Position) 
    /// computation. The result must depend only on time and positions (of any
    /// body or mobilizer), or earlier-stage state variables; it must not depend
    /// on any velocities. This Motion's own realizePosition() method will 
    /// already have been called. This will not be called if the u's are known 
    /// to be zero.
    virtual void calcPrescribedVelocity(const State& s, int nu, Real* u) const;

    /// Calculate the time derivative of the prescribed velocity. The udots 
    /// calculated here must be the exact time derivatives of the u's returned 
    /// by calcPrescribedVelocity(). So the calculation must be limited to the 
    /// same dependencies, plus the current value of this mobilizer's u's (or 
    /// the cross-mobilizer velocity V_FM because that depends only on those u's).
    ///
    /// This operator is called during the MatterSubsystem's realize(Dynamics) 
    /// computation. This Motion's own realizeDynamics() method will already 
    /// have been called. This will not be called if the udots are known to be 
    /// zero.
    virtual void calcPrescribedVelocityDot(const State& s, int nu, Real* udot) const;
    //@}

    /// @name Acceleration-only prescribed motion virtual
    ///
    /// This must be defined if the motion method is "Prescribed" and the 
    /// motion level is "Acceleration". In that case udot=udot(t,q,u) is 
    /// required. Note that Simbody passes in the number of u's (same as 
    /// number of udots) being prescribed; make sure you are seeing what you 
    /// expect.
    //@{
    /// This operator is called during the MatterSubsystem's realize(Dynamics) 
    /// computation. The result can depend on time, any positions, and any 
    /// velocities but must not depend on accelerations or reaction forces. This 
    /// Motion's own realizeDynamics() method will already have been called. 
    /// This will not be called if the udots are known to be zero.
    virtual void calcPrescribedAcceleration(const State& s, int nu, Real* udot) const;
    //@}

    /** 
     * @name Optional realize() virtual methods
     *
     * The following methods may optionally be overridden to do specialized 
     * realization for a Motion. These are called during the corresponding
     * realization stage of the containing MatterSubsystem.
     */
    //@{
    virtual void realizeTopology    (State&       state) const {}
    virtual void realizeModel       (State&       state) const {}
    virtual void realizeInstance    (const State& state) const {}
    virtual void realizeTime        (const State& state) const {}
    virtual void realizePosition    (const State& state) const {}
    virtual void realizeVelocity    (const State& state) const {}
    virtual void realizeDynamics    (const State& state) const {}
    virtual void realizeAcceleration(const State& state) const {}
    virtual void realizeReport      (const State& state) const {}
    //@}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOTION_H_
