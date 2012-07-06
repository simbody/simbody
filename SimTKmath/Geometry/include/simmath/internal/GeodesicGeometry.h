#ifndef SimTK_SIMBODY_GEODESICGEOMETRY_H_
#define SimTK_SIMBODY_GEODESICGEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKmath                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Ian Stavness, Michael Sherman                                     *
 * Contributors:                                                              *
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

/**
 * This is a temporary class for the geodesic specific
 * additions to ContactGeometry
 **/

#include "SimTKcommon.h"
#include "simmath/Differentiator.h"
#include "simmath/RungeKutta3Integrator.h"
#include "simmath/TimeStepper.h"
#include "simmath/internal/Geodesic.h"
#include "simmath/internal/ParticleOnSurfaceSystem.h"
#include "simmath/internal/BicubicSurface.h" // XXX compiler needed this
#include "simmath/internal/ContactGeometry.h"
#include <cassert>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}


namespace SimTK {


//XXX these magic constants need to be integrated into the class properly...
const Real startTime = 0;
const Real finalTime = Pi;
const Real integratorAccuracy = 1e-6;
const Real gamma = 0.01; // Baumgarte stabilization; params from Ascher1993
const Real alpha = gamma*gamma; // position stabilization
const Real beta = 2*gamma; // velocity stabilization


class SplitGeodesicError;

/*
 * A simple plane class
 */
class Plane {
public:
    Plane() : m_normal(1,0,0), m_offset(0) { }
    Plane(const Vec3& normal, const Real& offset) 
    :   m_normal(normal), m_offset(offset) { }

    Real getDistance(const Vec3& pt) const {
        return ~m_normal*pt - m_offset;
    }

    Vec3 getNormal() const {
        return m_normal;
    }

    Real getOffset() const {
        return m_offset;
    }

private:
    Vec3 m_normal;
    Real m_offset;
}; // class Plane


/**
 * A event handler to terminate integration when geodesic hits the plane.
 * For use with a ParticleAlongSurfaceSystem
 *
 **/
class GeodHitPlaneEvent : public TriggeredEventHandler {
public:
    GeodHitPlaneEvent()
    :   TriggeredEventHandler(Stage::Position) { }

    explicit GeodHitPlaneEvent(const Plane& plane)
    :   TriggeredEventHandler(Stage::Position) {
        m_plane = plane;
    }

    // event is triggered if distance of geodesic endpoint to plane is zero
    Real getValue(const State& state) const {
        Vec3 endpt(&state.getQ()[0]);
        Real dist =  m_plane.getDistance(endpt);
//        std::cout << "dist = " << dist << std::endl;
        return dist;
    }

    // This method is called whenever this event occurs.
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {

        // This should be triggered when geodesic endpoint to plane is zero.
        Vec3 endpt;
        const Vector& q = state.getQ();
        endpt[0] = q[0]; endpt[1] = q[1]; endpt[2] = q[2];
        Real dist = m_plane.getDistance(endpt);

        ASSERT(std::abs(dist) < 0.01 );
        shouldTerminate = true;
//        std::cout << "hit plane!" << std::endl;
    }

    void setPlane(const Plane& plane) const {
        m_plane = plane;
    }

    const Plane& getPlane() const {
        return m_plane;
    }

private:
    mutable Plane m_plane;

}; // class GeodHitPlaneEvent



class SimTK_SIMMATH_EXPORT GeodesicGeometry {

public:
    GeodesicGeometry(const ContactGeometry& geom) :
            geom(geom), ptOnSurfSys(geom, alpha, beta), splitGeodErr(0) {
        geodHitPlaneEvent = new GeodHitPlaneEvent();
        ptOnSurfSys.addEventHandler(geodHitPlaneEvent); // takes ownership
        ptOnSurfSys.realizeTopology();

    }


    /** Given two points, find a geodesic curve connecting them.
    If a preferred starting point is provided, find the geodesic curve that
    is closest to that point. Otherwise, find the shortest length geodesic.

    @param[in] xP            Coordinates of the first point.
    @param[in] xQ            Coordinates of the second point.
    @param[in] xSP           (Optional) Coordinates of a preferred point for the geodesic to be near
    @param[in] options       Parameters related to geodesic calculation
    @param[out] geod         On exit, this contains a geodesic between P and Q.
    **/
    void initGeodesic(const Vec3& xP, const Vec3& xQ, const Vec3& xSP,
            const GeodesicOptions& options, Geodesic& geod) const;


    /** Given two points and previous geodesic curve close to the points, find
    a geodesic curve connecting the points that is close to the previous geodesic.
    XXX if xP and xQ are the exact end-points of prevGeod; then geod = prevGeod;

    @param[in] xP            Coordinates of the first point.
    @param[in] xQ            Coordinates of the second point.
    @param[in] prevGeod      A previous geodesic that should be near the new one.
    @param[in] options       Parameters related to geodesic calculation
    @param[out] geod         On exit, this contains a geodesic between P and Q.
    **/
    void continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
            const GeodesicOptions& options, Geodesic& geod);


    /** Compute a geodesic curve of the given length, starting at the given point and
     in the given direction.

    @param[in] xP            Coordinates of the starting point for the geodesic.
    @param[in] tP            The starting tangent direction for the geodesic.
    @param[in] terminatingLength   The length that the resulting geodesic should have.
    @param[in] options       Parameters related to geodesic calculation
    @param[out] geod         On exit, this contains the calculated geodesic
    **/
    void calcGeodesicInDirectionUntilLengthReached(const Vec3& xP, const Vec3& tP,
            const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const;
    // XXX what to do if tP is not in the tangent plane at P -- project it?


    /** Compute a geodesic curve starting at the given point, starting in the given
     direction, and terminating at the given plane.

    @param[in] xP            Coordinates of the starting point for the geodesic.
    @param[in] tP            The starting tangent direction for the geodesic.
    @param[in] terminatingPlane   The plane in which the end point of the resulting geodesic should lie.
    @param[in] options       Parameters related to geodesic calculation
    @param[out] geod         On exit, this contains the calculated geodesic
    **/
    // XXX what to do if tP is not in the tangent plane at P -- project it?
    // XXX what to do if we don't hit the plane
   void calcGeodesicInDirectionUntilPlaneHit(const Vec3& P, const Vec3& tP,
            const Plane& terminatingPlane, const GeodesicOptions& options,
            Geodesic& geod) const;


    /** Utility method to find geodesic between P and Q with initial shooting
     directions tPhint and tQhint
     **/
    void calcGeodesic(const Vec3& xP, const Vec3& xQ,
            const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const;



    /**
     * Utility method to calculate the "geodesic error" between one geodesic
     * shot from P in the direction tP and another geodesic shot from Q in the
     * direction tQ
     **/
    static Vec2 calcGeodError(const GeodesicGeometry& gg,
            const Vec3& P, const Vec3& Q, const Vec3& tP, const Vec3& tQ);


    /** Temporary utility method for creating a continuous geodesic object
     *  from the two "split geodesics"
     *  XXX to be replaced by actual copy / append functionalily in Geodesic.h
     **/
    static void mergeGeodesics(const Geodesic& geodP, const Geodesic& geodQ,
            Geodesic& geod) {
        mergeGeodLists(geodP.getPoints(), geodQ.getPoints(), geod.updPoints());
        mergeGeodLists(geodP.getTangents(), geodQ.getTangents(), geod.updTangents());
        // TODO also merge other data, including arcLengths, etc.
    }

    static void mergeGeodLists(const Array_<Vec3>& geodP,
        const Array_<Vec3>& geodQ, Array_<Vec3>& geod) {
//        std::cout << "merge lists" << std::endl;

        int sizeP = geodP.size();
        int sizeQ = geodQ.size();
        geod.reserve(sizeP+sizeQ-1);

        for (int i = 0; i < sizeP-1; ++i) {
            geod.push_back(Vec3(geodP[i]));
        }
        geod.push_back(Vec3((geodP[sizeP-1]+geodQ[sizeQ-1])/2)); //midpt
        for (int i = sizeQ-2; i >= 0; --i) {
            geod.push_back(Vec3(geodQ[i]));
        }
    }


    /** Get the plane associated with the
        geodesic hit plane event handler  **/
    const Plane& getPlane() const {
        return geodHitPlaneEvent->getPlane();
    }

    /** Get the ContactGeometry object associated with this
        GeodesicGeometry object **/
    const ContactGeometry& getGeom() const {
        return geom;
    }

    /** @name Advanced/Obscure/Debugging **/
    /**@{**/

    /** Get the geodesic for access by visualizer **/
    const Geodesic& getGeodP() {
        return geodP;
    }

    /** Get the geodesic for access by visualizer **/
    const Geodesic& getGeodQ() {
        return geodQ;
    }
    /**@}**/

private:
    ContactGeometry geom;
    ParticleOnSurfaceSystem ptOnSurfSys;
    GeodHitPlaneEvent* geodHitPlaneEvent; // don't delete this

    class SplitGeodesicError; // local class
    mutable SplitGeodesicError* splitGeodErr;

    // temporary objects
    mutable Geodesic geodQ;
    mutable Geodesic geodP;

    friend class SplitGeodesicError;
}; // class GeodesicGeometry


/**
 * This class generates decoration for contact points and straight line path segments
 */
class PathDecorator : public DecorationGenerator {
public:
    PathDecorator(const Vector& x, const Vec3& O, const Vec3& I, const Vec3& color) :
            m_x(x), m_O(O), m_I(I), m_color(color) { }

    virtual void generateDecorations(const State& state,
            Array_<DecorativeGeometry>& geometry) {
//        m_system.realize(state, Stage::Position);

        Vec3 P, Q;
        P[0] = m_x[0]; P[1] = m_x[1]; P[2] = m_x[2];
        Q[0] = m_x[3]; Q[1] = m_x[4]; Q[2] = m_x[5];

        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(m_O));
        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(P));
        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(Q));
        geometry.push_back(DecorativeSphere(0.05).setColor(Black).setTransform(m_I));

        geometry.push_back(DecorativeLine(m_O,P)
                .setColor(m_color)
                .setLineThickness(2));
        geometry.push_back(DecorativeLine(Q,m_I)
                .setColor(m_color)
                .setLineThickness(2));

    }

private:
    const Vector& m_x; // x = ~[P Q]
    const Vec3& m_O;
    const Vec3& m_I;
    const Vec3& m_color;
    Rotation R_plane;
    Vec3 offset;
}; // class DecorationGenerator



/**
 * This class generates decoration for a plane
 */
class PlaneDecorator : public DecorationGenerator {
public:
    PlaneDecorator(const Plane& plane, const Vec3& color) :
            m_plane(plane), m_color(color) { }

    virtual void generateDecorations(const State& state,
            Array_<DecorativeGeometry>& geometry) {
//        m_system.realize(state, Stage::Position);

        // draw plane
        R_plane.setRotationFromOneAxis(UnitVec3(m_plane.getNormal()),
                CoordinateAxis::XCoordinateAxis());
        offset = 0;
        offset[0] = m_plane.getOffset();
        geometry.push_back(
                DecorativeBrick(Vec3(0.01,1,1))
                .setTransform(Transform(R_plane, R_plane*offset))
                .setColor(m_color)
                .setOpacity(.2));
    }

private:
    const Plane& m_plane;
    const Vec3& m_color;
    Rotation R_plane;
    Vec3 offset;
}; // class DecorationGenerator

} // namespace SimTK

#endif /*SimTK_SIMBODY_GEODESICGEOMETRY_H_*/
