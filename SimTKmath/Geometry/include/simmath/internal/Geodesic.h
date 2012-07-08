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
/**
 * This class stores quantities related to geodesic curves
 */

#include "SimTKcommon.h"

namespace SimTK {

class Geodesic {
public:
    /** Construct an empty geodesic. **/
    Geodesic() {clear();}

    Array_<Vec3>& updPoints() {
        return points;
    }

    const Array_<Vec3>& getPoints() const {
        return points;
    }

    void addPoint(const Vec3& P) {
        points.push_back(P);
    }


    Array_<Vec3>& updTangents() {
        return tangents;
    }

    const Array_<Vec3>& getTangents() const {
        return tangents;
    }

    void addTangent(const Vec3& tP) {
        tangents.push_back(tP);
    }

    Array_<Real>& updArcLengths() {
        return arcLengths;
    }

    const Array_<Real>& getArcLengths() const {
        return arcLengths;
    }

    void addArcLength(Real s) {
        arcLengths.push_back(s);
    }

    Real getLength() const {return arcLengths.empty() ? 0 : arcLengths.back();}

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
        points.clear(); 
        tangents.clear();
        arcLengths.clear();
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

private:
    Array_<Vec3> points;        // point at each knot
    Array_<Vec3> tangents;      // tangent at that point
    Array_<Real> arcLengths;    // arc length coord corresponding to that point

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
        const Array_<Vec3>& pts = m_geod.getPoints();
        Vec3 prevPt;
        for (int i = 0; i < (int) pts.size(); ++i) {
            Vec3 cur = pts[i];
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

#endif /*SimTK_SIMMATH_GEODESIC_H_*/
