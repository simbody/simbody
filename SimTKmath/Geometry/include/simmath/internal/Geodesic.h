#ifndef SimTK_SIMMATH_GEODESIC_H_
#define SimTK_SIMMATH_GEODESIC_H_

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

/** @file
This file defines the Geodesic class. **/

//==============================================================================
//                           GEODESIC CLASS
//==============================================================================


#include "SimTKcommon.h"
#include "simmath/internal/common.h"

namespace SimTK {

/** This class stores a geodesic curve after it has been determined. The curve
is represented by a discrete set of Frenet frames along its arc length, with
each frame providing a point on the curve, the tangent along the curve
there, surface normal, and binormal. The number of points is determined by
the accuracy to which the geodesic was calculated and the complexity of
the surface. For analytical geodesics, we'll sample the curve to generate
enough points for visualization but the accuracy will be machine precision.
The first point is at arclength s=0, the last is at the
actual length of the geodesic. We call the first point P and the last
point Q, and the geodesic arc length increases from P to Q, with the
tangent always pointing in the direction of increasing arc length. **/
class SimTK_SIMMATH_EXPORT Geodesic {
public:
    /** Construct an empty geodesic. **/
    Geodesic() {clear();}

    int getNumPoints() const {return (int)frenetFrames.size(); }

    /** Frenet frame of geodesic at arc length s:
        - origin: point on geodesic at s
        - z axis: outward surface (unit) normal n at s
        - y axis: unit tangent t at s in direction of increasing arc length
        - x axis: unit binormal at s: b=t X n (x=y X z)

    Note that this convention is different from that of a curve, because
    our normal vector always points outwards from the surface. **/
    const Array_<Transform>& getFrenetFrames() const {return frenetFrames;}
    Array_<Transform>&       updFrenetFrames() {return frenetFrames;}
    void addFrenetFrame(const Transform& Kf) {frenetFrames.push_back(Kf);}

    Array_<Real>& updArcLengths() {return arcLengths;}
    const Array_<Real>& getArcLengths() const {return arcLengths;}
    void addArcLength(Real s) {arcLengths.push_back(s);}

    Array_<Real>& updCurvatures() {return curvature;}
    const Array_<Real>& getCurvatures() const {return curvature;}
    void addCurvature(Real kappa) {curvature.push_back(kappa);}

    Array_<Vec2>& updDirectionalSensitivityPtoQ()
    {   return directionalSensitivityPtoQ; }
    const Array_<Vec2>& getDirectionalSensitivityPtoQ() const
    {   return directionalSensitivityPtoQ; }
    void addDirectionalSensitivityPtoQ(const Vec2& jP) {
        directionalSensitivityPtoQ.push_back(jP);
    }

    Array_<Vec2>& updDirectionalSensitivityQtoP()
    {   return directionalSensitivityQtoP; }
    const Array_<Vec2>& getDirectionalSensitivityQtoP() const
    {   return directionalSensitivityQtoP; }
    void addDirectionalSensitivityQtoP(const Vec2& jQ) {
        directionalSensitivityQtoP.push_back(jQ);
    }

    Array_<Vec2>& updPositionalSensitivityPtoQ()
    {   return positionalSensitivityPtoQ; }
    const Array_<Vec2>& getPositionalSensitivityPtoQ() const
    {   return positionalSensitivityPtoQ; }
    void addPositionalSensitivityPtoQ(const Vec2& jtP) {
        positionalSensitivityPtoQ.push_back(jtP);
    }

    Array_<Vec2>& updPositionalSensitivityQtoP()
    {   return positionalSensitivityQtoP; }
    const Array_<Vec2>& getPositionalSensitivityQtoP() const
    {   return positionalSensitivityQtoP; }
    void addPositionalSensitivityQtoP(const Vec2& jtQ) {
        positionalSensitivityQtoP.push_back(jtQ);
    }

    void setTorsionAtP(Real tauP) {torsionAtP = tauP;}
    void setTorsionAtQ(Real tauQ) {torsionAtQ = tauQ;}
    void setBinormalCurvatureAtP(Real muP) {binormalCurvatureAtP = muP;}
    void setBinormalCurvatureAtQ(Real muQ) {binormalCurvatureAtQ = muQ;}

    /** Return the total arc length of this geodesic curve. Will return zero
    if no curve has been calculated. **/
    Real getLength() const {return arcLengths.empty() ? 0 : arcLengths.back();}

    /** Given the time derivatives of the surface coordinates of P and Q,
    calculate the rate of change of length of this geodesic. Only the
    components of \a xdotP and \a xdotQ along the curve tangent directions
    can affect its length; moves in the normal or binormal direction leave
    the length unchanged. **/
    Real calcLengthDot(const Vec3& xdotP, const Vec3& xdotQ) const
    {   return ~xdotQ*getTangentQ() - ~xdotP*getTangentP(); }

    /** Return the location on the surface of the geodesic's starting point P,
    measured and expressed in the surface frame S. **/
    const Vec3& getPointP() const {return frenetFrames.front().p();}
    /** Return the location on the surface of the geodesic's ending point Q,
    measured and expressed in the surface frame S. **/
    const Vec3& getPointQ() const {return frenetFrames.back().p();}

    /** Return the surface outward unit normal at P, which is aligned
    with the curve normal there but will have opposite sign if the
    geodesic curvature is positive at P. **/
    const UnitVec3& getNormalP() const {return frenetFrames.front().z();}
    /** Return the surface outward unit normal at Q, which is aligned
    with the curve normal there but will have opposite sign if the
    geodesic curvature is positive at Q. **/
    const UnitVec3& getNormalQ() const {return frenetFrames.back().z();}

    /** Return the unit tangent to the geodesic at P, pointing in the
    direction of increasing arc length parameters (i.e., towards Q). **/
    const UnitVec3& getTangentP() const {return frenetFrames.front().y();}
    /** Return the unit tangent to the geodesic at Q, pointing in the
    direction of increasing arc length parameters (i.e., away from P). **/
    const UnitVec3& getTangentQ() const {return frenetFrames.back().y();}

    /** Return the unit binormal vector to the curve at P, defined as
    bP = tP X nP. **/
    const UnitVec3& getBinormalP() const {return frenetFrames.front().x();}
    /** Return the unit binormal vector to the curve at Q, defined as
    bQ = tQ X nQ. **/
    const UnitVec3& getBinormalQ() const {return frenetFrames.back().x();}

    /** Return the geodesic normal curvature at P, defined to be positive when
    the surface is convex in the curve tangent direction at P, negative if the
    surface is concave in that direction. This is a scalar
    kappaP=-dot(DtP/ds,nP). Note that the geodesic curvature is
    the same as the surface curvature in the curve tangent direction. **/
    Real getCurvatureP() const {return curvature.front();}
    /** Return the geodesic normal curvature at Q, defined to be positive when
    the surface is convex in the curve tangent direction at Q, negative if the
    surface is concave in that direction. This is a scalar
    kappaQ=-dot(DtQ/ds,nQ). Note that the geodesic curvature is
    the same as the surface curvature in the curve tangent direction, and
    remember that the curve tangent at Q points \e away from P, that is, off
    the end of the geodesic. **/
    Real getCurvatureQ() const {return curvature.back();}

    /** Return the geodesic torsion at P, that is, the twisting of the
    Frenet frame as you move along the tangent towards Q. The sign follows
    the right hand rule with your thumb directed along the tangent. This
    is a scalar tauP=-dot(DbP/Ds,nP). **/
    Real getTorsionP() const {return torsionAtP;}
    /** Return the geodesic torsion at Q, that is, the twisting of the
    Frenet frame as you move along the tangent away from P, that is, off the
    end of the geodesic. The sign follows the right hand rule with your thumb
    directed along the tangent. This is a scalar tauQ=-dot(DbQ/Ds,nQ). **/
    Real getTorsionQ() const {return torsionAtQ;}

    /** Return the \e surface curvature in the binormal direction at P; don't
    confuse this with the geodesic torsion at P. Surface curvature in a
    direction is a property of the surface independent of any curve. **/
    Real getBinormalCurvatureP() const {return binormalCurvatureAtP;}
    /** Return the \e surface curvature in the binormal direction at Q; don't
    confuse this with the geodesic torsion at Q. Surface curvature in a
    direction is a property of the surface independent of any curve. **/
    Real getBinormalCurvatureQ() const {return binormalCurvatureAtQ;}

    /** Return jP, the Jacobi field term giving the sensitivity of the P end
    of the geodesic to changes in tangent direction at the Q end, assuming
    the geodesic length is fixed. This is a scalar jP=dot(DP/DthetaQ,bP) where
    thetaQ is a right-hand-rule rotation about the normal nQ at Q. Caution:
    jP and jQ have opposite signs. **/
    Real getJacobiP() const {return directionalSensitivityQtoP.front()[0];}
    /** Return jQ, the Jacobi field term giving the sensitivity of the Q end
    of the geodesic to changes in tangent direction at the P end, assuming
    the geodesic length is fixed. This is a scalar jQ=-dot(DQ/DthetaP,bQ)
    where thetaP is a right-hand-rule rotation about the normal nP at P.
    Note that jQ has the opposite sign from jP. **/
    Real getJacobiQ() const {return directionalSensitivityPtoQ.back()[0];}

    /** Return the derivative of jP with respect to s, the arc length of the
    geodesic (which always runs from P to Q). That is jPDot = d/ds jP. **/
    // Note sign change here -- we calculated this by integrating backwards
    // so the arc length we used had the opposite sign from the true arc
    // length parameter.
    Real getJacobiPDot() const {return -directionalSensitivityQtoP.front()[1];}
    /** Return the derivative of jQ with respect to s, the arc length of the
    geodesic. That is, jQDot = d/ds jQ. **/
    Real getJacobiQDot() const {return directionalSensitivityPtoQ.back()[1];}

    // XXX testing
    Real getJacobiTransP() const {return positionalSensitivityQtoP.front()[0];}
    Real getJacobiTransQ() const {return positionalSensitivityPtoQ.back()[0];}
    Real getJacobiTransPDot() const {return -positionalSensitivityQtoP.front()[1];}
    Real getJacobiTransQDot() const {return positionalSensitivityPtoQ.back()[1];}

    /** Clear the data in this geodesic, returning it to its default-constructed
    state, although memory remains allocated. **/
    void clear() {
        arcLengths.clear();
        frenetFrames.clear();
        directionalSensitivityPtoQ.clear();
        directionalSensitivityQtoP.clear();
        positionalSensitivityPtoQ.clear();
        positionalSensitivityQtoP.clear();
        curvature.clear();
        torsionAtP = torsionAtQ = NaN;
        binormalCurvatureAtP = binormalCurvatureAtQ = NaN;
        convexFlag = shortestFlag = false;
        initialStepSizeHint = achievedAccuracy = NaN;
    }

    void setIsConvex(bool isConvex) {convexFlag = isConvex;}
    void setIsShortest(bool isShortest) {shortestFlag = isShortest;}
    void setInitialStepSizeHint(Real sz) {initialStepSizeHint=sz;}
    void setAchievedAccuracy(Real acc) {achievedAccuracy=acc;}

    bool isConvex() const {return convexFlag;}
    bool isShortest() const {return shortestFlag;}
    Real getInitialStepSizeHint() const {return initialStepSizeHint;}
    Real getAchievedAccuracy() const {return achievedAccuracy;}

    void dump(std::ostream& o) const;

private:
    // All these arrays are the same length when the geodesic is complete.
    Array_<Real>      arcLengths; // arc length coord corresponding to point
    Array_<Transform> frenetFrames; // see above for more info
    Array_<Vec2>      directionalSensitivityPtoQ; // jQ and jQdot
    Array_<Vec2>      directionalSensitivityQtoP; // jP and -jPdot
    Array_<Vec2>      positionalSensitivityPtoQ; // jtQ and jtQdot
    Array_<Vec2>      positionalSensitivityQtoP; // jtP and -jtPdot
    Array_<Real>      curvature; // normal curvature kappa in tangent direction
    // These are only calculated at the end points.
    Real              torsionAtP, torsionAtQ; // torsion tau (only at ends)
    Real              binormalCurvatureAtP, binormalCurvatureAtQ;

    // This flag is set true if curvature[i]>=0 for all i.
    bool convexFlag; // is this geodesic over a convex surface?

    bool shortestFlag; // XXX is this geodesic the shortest one of the surface?
    Real initialStepSizeHint; // the initial step size to be tried when integrating this geodesic
    Real achievedAccuracy; // the accuracy to which this geodesic curve has been calculated
};


/**
 * This class generates decoration (line segments) for a geodesic curve
 */
class GeodesicDecorator : public DecorationGenerator {
public:
    GeodesicDecorator(const Geodesic& geod, const Vec3& color) :
            m_geod(geod), m_color(color) { }

    virtual void generateDecorations(const State& state,
            Array_<DecorativeGeometry>& geometry) {
//        m_system.realize(state, Stage::Position);

        // draw connected line segments from pts
        const Array_<Transform>& Kfs = m_geod.getFrenetFrames();
        Vec3 prevPt;
        for (int i = 0; i < (int) Kfs.size(); ++i) {
            Vec3 cur = Kfs[i].p();
            if (i > 0) {
                geometry.push_back(
                        DecorativeLine(prevPt, cur)
                        .setColor(m_color)
                        .setLineThickness(3));
            }
            prevPt = cur;

            geometry.push_back(DecorativeFrame(Real(.2)).setTransform(Kfs[i]));
        }
    }

private:
    const Geodesic& m_geod;
    const Vec3& m_color;
};




/**
 * This class stores options for calculating geodesics
 */
class GeodesicOptions {
    // XXX stub
};


} // namespace SimTK

#endif // SimTK_SIMMATH_GEODESIC_H_
