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
        - x axis: unit tangent t at s in direction of increasing arc length
        - y axis: unit binormal at s: b=n X t (y=z X x) 
    **/
    const Array_<Transform>& getFrenetFrames() const {return frenetFrames;}
    Array_<Transform>&       updFrenetFrames() {return frenetFrames;}

    void addFrenetFrame(const Transform& Kf) {frenetFrames.push_back(Kf);}

    Array_<Real>& updArcLengths() {return arcLengths;}
    const Array_<Real>& getArcLengths() const {return arcLengths;}

    void addArcLength(Real s) {arcLengths.push_back(s);}

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

    /** Return the total arc length of this geodesic curve. Will return zero
    if no curve has been calculated. **/
    Real getLength() const {return arcLengths.empty() ? 0 : arcLengths.back();}

    const Vec3& getPointP() const {return frenetFrames.front().p();}
    const Vec3& getPointQ() const {return frenetFrames.back().p();}
    const UnitVec3& getNormalP() const {return frenetFrames.front().z();}
    const UnitVec3& getNormalQ() const {return frenetFrames.back().z();}
    const UnitVec3& getTangentP() const {return frenetFrames.front().x();}
    const UnitVec3& getTangentQ() const {return frenetFrames.back().x();}
    const UnitVec3& getBinormalP() const {return frenetFrames.front().y();}
    const UnitVec3& getBinormalQ() const {return frenetFrames.back().y();}
    Real getSensitivityP() const {return directionalSensitivityPtoQ.back()[0];}
    Real getSensitivityQ() const {return directionalSensitivityQtoP.front()[0];}


    /** TODO: Given the time derivatives of the surface coordinates of P and Q,
    calculate the rate of change of length of this geodesic. **/
    Real calcLengthDot(const Vec3& xdotP, const Vec3& xdotQ) const 
    {   return 0; }

    // ARC LENGTH METHODS

    /** Given arc length coordinate return the corresponding geodesic point. **/
    Vec3 findPoint(Real s) const;
    Vec3 findTangent(Real s) const;
    Vec3 findNormal(Real s) const;
    Vec3 findBinormal(Real s) const; // tangent X normal

    /** Clear the data in this geodesic, returning it to its default-constructed
    state, although memory remains allocated. **/
    void clear() {
        arcLengths.clear();
        frenetFrames.clear(); 
        directionalSensitivityPtoQ.clear(); 
        directionalSensitivityQtoP.clear(); 
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
    Array_<Real>      arcLengths; // arc length coord corresponding to that point
    Array_<Transform> frenetFrames;
    Array_<Vec2>      directionalSensitivityPtoQ; // j and jdot
    Array_<Vec2>      directionalSensitivityQtoP;


    // XXX other members:
    bool convexFlag; // is this geodesic over a convex surface?
    bool shortestFlag; // is this geodesic the shortest one of the surface?
    Real initialStepSizeHint; // the initial step size to be tried when integrating this geodesic
    Real achievedAccuracy; // the accuracy to which this geodesic curve has been calculated
//    Array_<JacobiFieldInfo> partials;
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
