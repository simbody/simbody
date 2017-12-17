#ifndef SimTK_SIMMATH_CONTACT_GEOMETRY_IMPL_H_
#define SimTK_SIMMATH_CONTACT_GEOMETRY_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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


#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/OBBTree.h"
#include "simmath/internal/ParticleConSurfaceSystem.h"
#include "simmath/Differentiator.h"
#include "simmath/internal/ContactGeometry.h"

#include <atomic>
#include <limits>

namespace SimTK {

class SplitGeodesicError;


//==============================================================================
//                             CONTACT GEOMETRY IMPL
//==============================================================================
class SimTK_SIMMATH_EXPORT ContactGeometryImpl {
public:
    ContactGeometryImpl() 
    :   myHandle(0), ptOnSurfSys(0), geodHitPlaneEvent(0), vizReporter(0),
        splitGeodErr(0), numGeodesicsShot(0)
    {
        createParticleOnSurfaceSystem();
    }
    ContactGeometryImpl(const ContactGeometryImpl& source)
    :   myHandle(0), ptOnSurfSys(0), geodHitPlaneEvent(0), vizReporter(0),
        splitGeodErr(0), numGeodesicsShot(0) 
    {}

    virtual ~ContactGeometryImpl() {
        clearMyHandle();
        clearParticleOnSurfaceSystem();
    }

    /* Create a new ContactGeometryTypeId and return this unique integer 
    (thread safe). Each distinct type of ContactGeometry should use this to
    initialize a static variable for that concrete class. */
    static ContactGeometryTypeId  createNewContactGeometryTypeId()
    {   static std::atomic<int> nextAvailableId(1);
        return ContactGeometryTypeId(nextAvailableId++); }

    virtual ContactGeometryImpl*  clone() const = 0;
    virtual ContactGeometryTypeId getTypeId() const = 0;

    virtual DecorativeGeometry createDecorativeGeometry() const {
        SimTK_THROW2(Exception::UnimplementedVirtualMethod,
                "ContactGeometryImpl", "createDecorativeGeometry()");
    }

    virtual Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                                  UnitVec3& normal) const = 0;

    virtual bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                               Real& distance, UnitVec3& normal) const = 0;

    virtual void getBoundingSphere(Vec3& center, Real& radius) const = 0;


    virtual bool isSmooth() const = 0;
    virtual bool isConvex() const = 0;
    virtual bool isFinite() const = 0;

    // Smooth surfaces only.
    virtual void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "ContactGeometryImpl", "calcCurvature"); }

    // Smooth surfaces only.
    virtual const Function& getImplicitFunction() const
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "ContactGeometryImpl", "getImplicitFunction"); }

    // Convex surfaces only.
    virtual Vec3 calcSupportPoint(const UnitVec3& direction) const
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "ContactGeometryImpl", "calcSupportPoint"); }

    const OBBTree& getOBBTree() const {return obbTree;}


    Real  calcSurfaceValue(const Vec3& point) const;
    UnitVec3 calcSurfaceUnitNormal(const Vec3& point) const;
    Vec3  calcSurfaceGradient(const Vec3& point) const;
    Mat33 calcSurfaceHessian(const Vec3& point) const;
    Real  calcGaussianCurvature(const Vec3& gradient,
                                const Mat33& Hessian) const;
    Real  calcSurfaceCurvatureInDirection(const Vec3& point, 
                                          const UnitVec3& direction) const;
    // Generic method for calculating principal curvatures kmax,kmin and
    // corresponding unit tangent vector directions R_SP.x() and R_SP.y().
    // R_SP.z() is the surface unit normal at P, with z=x X y.
    void calcSurfacePrincipalCurvatures(const Vec3& point,
                                        Vec2& curvature,
                                        Rotation& R_SP) const;

    Vec3 projectDownhillToNearestPoint(const Vec3& Q) const;

    bool trackSeparationFromLine(const Vec3& pointOnLine,
                        const UnitVec3& directionOfLine,
                        const Vec3& startingGuessForClosestPoint,
                        Vec3& newClosestPointOnSurface,
                        Vec3& closestPointOnLine,
                        Real& height) const;

    // Geodesic evaluators


    // Given two points, find a geodesic curve connecting them.
    void initGeodesic(const Vec3& xP, const Vec3& xQ, const Vec3& xSP,
            const GeodesicOptions& options, Geodesic& geod) const;


    // Given two points and previous geodesic curve close to the points, find
    // a geodesic curve connecting the points that is close to the previous geodesic.
    void continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
            const GeodesicOptions& options, Geodesic& geod) const;

    // Given two points (which should be close together) create a two-point
    // geodesic that is a straight line between the points.
    void makeStraightLineGeodesic(const Vec3& xP, const Vec3& xQ,
            const UnitVec3& defaultDirectionIfNeeded,
            const GeodesicOptions& options, Geodesic& geod) const;

    // Compute a geodesic curve starting at the given point, starting in the
    // given direction, and terminating at the given length.
    void shootGeodesicInDirectionUntilLengthReached(const Vec3& xP, const UnitVec3& tP,
            const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const;

    // Compute a geodesic curve starting at the given point, starting in the
    // given direction, and terminating when it hits the given plane.
    void shootGeodesicInDirectionUntilPlaneHit(const Vec3& xP, const UnitVec3& tP,
            const Plane& terminatingPlane, const GeodesicOptions& options,
            Geodesic& geod) const;

    // Utility method to find geodesic between P and Q with initial shooting
    // directions tPhint and tQhint
    void calcGeodesic(const Vec3& xP, const Vec3& xQ,
            const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const;
    
    // Use a method that generates orthogonal error conditions in the tangent
    // and binormal directions.
    void calcGeodesicUsingOrthogonalMethod
       (const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, Real lengthHint, Geodesic& geod) const; 

    // Utility method to calculate the "geodesic error" between one geodesic
    // shot from P in the direction tP and another geodesic shot from Q in the
    // direction tQ.
    Vec2 calcSplitGeodError(const Vec3& xP, const Vec3& xQ,
                       const UnitVec3& tP, const UnitVec3& tQ,
                       Geodesic* geod=0) const;

    // analytical versions of the geodesic API methods

    virtual void shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
            const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const;

    virtual void shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
            const Plane& terminatingPlane, const GeodesicOptions& options,
            Geodesic& geod) const;

    virtual void calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
            const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const;

    Vec2 calcSplitGeodErrorAnalytical(const Vec3& P, const Vec3& Q,
                       const UnitVec3& tP, const UnitVec3& tQ,
                       Geodesic* geod=0) const;


    // Utility method to calculate the "geodesic error" between the end-points
    // of two geodesics.
    Vec2 calcError(const Geodesic& geodP, const Geodesic& geodQ) const;


    // Calculate a unit tangent vector based on angle theta measured from the
    // binormal axis of the given rotation matrix
    static UnitVec3 calcUnitTangentVec(const Real& theta, const Rotation& R_GS) {
        UnitVec3 tR_S(std::sin(theta), std::cos(theta), 0);
        return R_GS*tR_S;
    }

    // Temporary utility method for creating a continuous geodesic object
    //  from the two "split geodesics"
    //  XXX to be replaced by actual copy / append functionalily in Geodesic.h
    static void mergeGeodesics(const Geodesic& geodP, const Geodesic& geodQ,
            Geodesic& geod) {

        mergeGeodFrames(geodP.getFrenetFrames(), geodQ.getFrenetFrames(),
                        geod.updFrenetFrames());
        mergeArcLengths(geodP.getArcLengths(), geodQ.getArcLengths(),
                        geod.updArcLengths());
        mergeJacobiField(geodP.getDirectionalSensitivityPtoQ(), 
                         geodQ.getDirectionalSensitivityQtoP(),
                         geod.updDirectionalSensitivityPtoQ());
        mergeJacobiField(geodP.getDirectionalSensitivityQtoP(), 
                         geodQ.getDirectionalSensitivityPtoQ(),
                         geod.updDirectionalSensitivityQtoP());
        // TODO: calculate QtoP
        geod.setIsConvex(geodP.isConvex() && geodQ.isConvex());
        // TODO: what to do with "is shortest"?
        // Take the smaller of the two suggested step sizes.
        geod.setInitialStepSizeHint(std::min(geodP.getInitialStepSizeHint(),
                                             geodQ.getInitialStepSizeHint()));
        // Take the larger (sloppier) of the two accuracies.
        geod.setAchievedAccuracy(std::max(geodP.getAchievedAccuracy(),
                                          geodQ.getAchievedAccuracy()));

    }

    // Temporary utility method for joining geodP and the reverse of geodQ.
    // Optionally negate the Q vector (used for tangents). The last entries
    // for both geodesics are supposed to represent the same point on the
    // combined geodesic so we'll average them. The final geodesic will
    // thus have one fewer point than the sum of the two half-geodesics.
    // Note convention for Transform<->Frenet Frame mapping:
    //    z = normal, y = tangent, x = binormal (tXn)
    static void mergeGeodFrames(const Array_<Transform>& geodP,
                                const Array_<Transform>& geodQ,
                                Array_<Transform>&       geod)
    {
        const int sizeP = geodP.size(), sizeQ = geodQ.size();
        assert(sizeP > 0 && sizeQ > 0);
        geod.clear();
        geod.reserve(sizeP+sizeQ-1);
        geod = geodP; // the first part is the same as P

        const Transform& Pf = geodP.back();
        const Transform& Qf = geodQ.back();
        // Recalculate the midpoint by averaging the ends.
        Vec3     midpoint = (Pf.p() + Qf.p())/2;
        // TODO: get exact normal at midpoint from surface
        UnitVec3 midnormal = UnitVec3((Pf.z() + Qf.z())/2); 
        Vec3     midtangent = (Pf.y() - Qf.y())/2; // approx

        // Replace the midpoint Frenet frame with one at the average of
        // the P and Q endpoints, using the calculated normal as the z axis,
        // and the average tangent as y axis. We'll let Rotation 
        // perpendicularize this, so the actual y axis may deviate slightly
        // from the calculated average normal. Then b=y X z is the x axis.
        geod.back() = Transform(Rotation(midnormal, ZAxis, midtangent, YAxis),
                                midpoint);

        // Now run through the Q-side geodesic backwards reversing its
        // tangent (y) and binormal (x) vectors while preserving the point
        // and surface normal. That's a 180 degree rotation about z so it
        // is still right handed!
        for (int i = sizeQ-2; i >= 0; --i) {
            const Transform& F = geodQ[i]; // old Frenet frame
            Transform Frev; // Same, but with x,y directions negated
            Frev.updR().setRotationFromUnitVecsTrustMe(-F.x(), -F.y(), F.z());
            Frev.updP() = F.p();
            geod.push_back(Frev);
        }
    }

    // Temporary utility for merging the arc lengths (knot coordinates) for
    // two half-geodesics. The last entry in each list is supposed to be the
    // same point.
    static void mergeArcLengths(const Array_<Real>& knotsP,
                                const Array_<Real>& knotsQ,
                                Array_<Real>&       knots)
    {
        const int sizeP = knotsP.size(), sizeQ = knotsQ.size();
        assert(sizeP > 0 && sizeQ > 0);
        knots.clear();
        knots.reserve(sizeP+sizeQ-1);
        knots = knotsP; // the first part is the same as P

        const Real lengthP = knotsP.back(), lengthQ = knotsQ.back();
        for (int i = sizeQ-2; i >= 0; --i)
            knots.push_back(lengthP + (lengthQ-knotsQ[i]));
    }

    static void mergeJacobiField(const Array_<Vec2>& jP,
                                 const Array_<Vec2>& jQ,
                                 Array_<Vec2>&       j)
    {
        const int sizeP = jP.size(), sizeQ = jQ.size();
        assert(sizeP > 0 && sizeQ > 0);
        j.clear();
        j.reserve(sizeP+sizeQ-1);
        j = jP; // the first part is the same as P
        //TODO: this is the wrong-direction field!
        for (int i = sizeQ-2; i >= 0; --i)
            j.push_back(jQ[i]);
    }


    // Utility method to used by calcGeodesicInDirectionUntilPlaneHit
    //  and calcGeodesicInDirectionUntilLengthReached
    void shootGeodesicInDirection(const Vec3& P, const UnitVec3& tP,
            const Real& finalTime, const GeodesicOptions& options,
            Geodesic& geod) const;
    void shootGeodesicInDirection2(const Vec3& P, const UnitVec3& tP,
            const Real& finalTime, const GeodesicOptions& options,
            Geodesic& geod) const;

    // Utility method to integrate geodesic backwards to fill in the Q to P
    // directional sensitivity.
    void calcGeodesicReverseSensitivity
       (Geodesic& geod, const Vec2& initJacobi) const;
    void calcGeodesicReverseSensitivity2
       (Geodesic& geod, const Vec2& initJRot, const Vec2& initJTrans) const;

    // Utility method to calculate the "geodesic error" between one geodesic
    // shot from P in the direction thetaP and another geodesic shot from Q in the
    // direction thetaQ given the pre-calculated member basis R_SP, R_SQ.
    // We optionally return the resulting "kinked" geodesic, which is the real
    // one if the returned errors are below tolerance.
    Vec2 calcSplitGeodError(const Vec3& xP, const Vec3& xQ,
                       Real thetaP, Real thetaQ,
                       Geodesic* geodesic=0) const;

    // Utility method to calculate the "orthogonal error" between a geodesic
    // shot from P in the direction thetaP with given length and a point Q.
    // We are actually shooting the geodesic here and we'll return it in
    // the indicated variable. Note that we *do not* do the backwards 
    // integration for the reverse directional sensitivity; if you want to
    // return this as the final geodesic be sure to fill in that term.
    Vec2 calcOrthogonalGeodError(const Vec3& xP, const Vec3& xQ,
                       Real thetaP, Real length,
                       Geodesic& geodesic) const;

    // Utility method to calculate the "geodesic error jacobian" between one geodesic
    // shot from P in the direction thetaP and another geodesic shot from Q in the
    // direction thetaQ given the pre-calculated member basis R_SP, R_SQ
    Mat22 calcSplitGeodErrorJacobian(const Vec3& P, const Vec3& Q,
            const Real& thetaP, const Real& thetaQ, Differentiator::Method method) const;

    // Compute rotation matrix using the normal at the given point and the
    // given direction.
    Rotation calcTangentBasis(const Vec3& point, const Vec3& dir) const {
        const UnitVec3 n = calcSurfaceUnitNormal(point);
        Rotation R_GS;
        R_GS.setRotationFromTwoAxes(n, ZAxis, dir, XAxis);
        return R_GS;
    }

    // Get the plane associated with the geodesic hit plane event handler
    const Plane& getPlane() const {
        return geodHitPlaneEvent->getPlane();
    }

    // Set the plane associated with the geodesic hit plane event handler
    void setPlane(const Plane& plane) const {
        return geodHitPlaneEvent->setPlane(plane);
    }


    // Get the geodesic for access by visualizer
    const Geodesic& getGeodP() const {
        return geodP;
    }

    // Get the geodesic for access by visualizer
    const Geodesic& getGeodQ() const {
        return geodQ;
    }

    const int getNumGeodesicsShot() const {
        return numGeodesicsShot;
    }

    void addVizReporter(ScheduledEventReporter* reporter) const {
        vizReporter = reporter;
        ptOnSurfSys->addEventReporter(vizReporter); // takes ownership
        ptOnSurfSys->realizeTopology();
    }

    class SplitGeodesicError; // local class
    friend class SplitGeodesicError;
    class OrthoGeodesicError; // local class
    friend class OrthoGeodesicError;

    void createParticleOnSurfaceSystem() {
        ptOnSurfSys = new ParticleConSurfaceSystem(*this);
        geodHitPlaneEvent = new GeodHitPlaneEvent();
        ptOnSurfSys->addEventHandler(geodHitPlaneEvent); // takes ownership
        ptOnSurfSys->realizeTopology();
    }

    void clearParticleOnSurfaceSystem() {
        delete ptOnSurfSys;
        ptOnSurfSys = 0;
        geodHitPlaneEvent = 0; // was deleted by the system
    }

    ContactGeometry* getMyHandle() {return myHandle;}
    void setMyHandle(ContactGeometry& h) {myHandle = &h;}
    void clearMyHandle() {myHandle = 0;}
protected:
    ContactGeometry*        myHandle;
    OBBTree                 obbTree;

    ParticleConSurfaceSystem* ptOnSurfSys;
    GeodHitPlaneEvent* geodHitPlaneEvent; // don't delete this
    mutable ScheduledEventReporter* vizReporter; // don't delete this
    mutable SplitGeodesicError* splitGeodErr;

    // temporary objects
    mutable Geodesic geodP;
    mutable Geodesic geodQ;
    mutable Rotation R_SP;
    mutable Rotation R_SQ;
    mutable int numGeodesicsShot;

};



//==============================================================================
//                             HALF SPACE IMPL
//==============================================================================
class HalfSpaceImplicitFunction : public Function {
public:
    HalfSpaceImplicitFunction() : ownerp(0) {}
    HalfSpaceImplicitFunction(const ContactGeometry::HalfSpace::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::HalfSpace::Impl& owner) {ownerp=&owner;}

    // Value is positive for x>0.
    Real calcValue(const Vector& P) const override {return P[0];}
    // First derivative w.r.t. x is 1, all else is zero.
    Real calcDerivative(const Array_<int>& components, 
                        const Vector& P) const override
    {   if (components.empty()) return calcValue(P);
        if (components.size()==1 && components[0]==0) return 1;
        return 0; }

    int getArgumentSize() const override {return 3;}
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::HalfSpace::Impl* ownerp; // just a ref.; don't delete
};


class ContactGeometry::HalfSpace::Impl : public ContactGeometryImpl {
public:
    Impl() : ContactGeometryImpl() {
    }
    ContactGeometryImpl* clone() const override {
        return new Impl();
    }

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const override;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    bool isSmooth() const override {return true;}
    bool isConvex() const override {return false;}
    bool isFinite() const override {return false;}

    // Curvature is zero everywhere. Since the half plane occupies x>0 in
    // its own frame, the surface normal is -x, and -x,y,-z forms a right 
    // handed set.
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const override
    {   curvature = 0;
        orientation.setRotationFromUnitVecsTrustMe
            (UnitVec3(-XAxis), UnitVec3(YAxis), UnitVec3(-ZAxis));
    }

    const Function& getImplicitFunction() const override {return function;}

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    HalfSpaceImplicitFunction function;
};



//==============================================================================
//                              CYLINDER IMPL
//==============================================================================
class CylinderImplicitFunction : public Function {
public:
    CylinderImplicitFunction() : ownerp(0) {}
    CylinderImplicitFunction(const ContactGeometry::Cylinder::Impl& owner)
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::Cylinder::Impl& owner) {ownerp=&owner;}
    Real calcValue(const Vector& x) const override;
    Real calcDerivative(const Array_<int>& derivComponents,
                        const Vector& x) const override;
    int getArgumentSize() const override {return 3;}
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::Cylinder::Impl* ownerp; // just a reference; don't delete
};

class ContactGeometry::Cylinder::Impl : public ContactGeometryImpl {
public:
    explicit Impl(Real radius) : radius(radius) {
        function.setOwner(*this);
    }

    ContactGeometryImpl* clone() const override {
        return new Impl(radius);
    }
    Real getRadius() const {
        return radius;
    }
    void setRadius(Real r) {
        radius = r;
    }

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside,
                          UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction,
                       Real& distance, UnitVec3& normal) const override;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    bool isSmooth() const override {return true;} //TODO: only for infinite
    bool isConvex() const override {return true;}
    bool isFinite() const override {return false;}

    Vec3 calcSupportPoint(const UnitVec3& direction) const override {
        assert(false);
        return Vec3(NaN);
    }

    void calcCurvature(const Vec3& point, Vec2& curvature,
                       Rotation& orientation) const override;

    void shootGeodesicInDirectionUntilLengthReachedAnalytical
       (const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, 
        Geodesic& geod) const override;

    void shootGeodesicInDirectionUntilPlaneHitAnalytical
       (const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const override;

    void calcGeodesicAnalytical
       (const Vec3& xP, const Vec3& xQ, const Vec3& tPhint, const Vec3& tQhint, 
        Geodesic& geod) const override;

    const Function& getImplicitFunction() const override {
        return function;
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id =
            createNewContactGeometryTypeId();
        return id;
    }
private:
    Real                    radius;
    CylinderImplicitFunction  function;
};



//==============================================================================
//                                SPHERE IMPL
//==============================================================================
class SphereImplicitFunction : public Function {
public:
    SphereImplicitFunction() : ownerp(0) {}
    SphereImplicitFunction(const ContactGeometry::Sphere::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::Sphere::Impl& owner) {ownerp=&owner;}
    Real calcValue(const Vector& x) const override;
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vector& x) const override;
    int getArgumentSize() const override {return 3;}
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::Sphere::Impl* ownerp; // just a reference; don't delete
};

class ContactGeometry::Sphere::Impl : public ContactGeometryImpl {
public:
    explicit Impl(Real radius) : radius(radius) {
        function.setOwner(*this);
        createOBBTree(); 
    }

    ContactGeometryImpl* clone() const override {
        return new Impl(radius);
    }
    Real getRadius() const {
        return radius;
    }
    void setRadius(Real r) {
        radius = r;
    }

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const override;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    bool isSmooth() const override {return true;}
    bool isConvex() const override {return true;}
    bool isFinite() const override {return true;}

    Vec3 calcSupportPoint(const UnitVec3& direction) const override {
        return radius*direction;
    }
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const override;

    void shootGeodesicInDirectionUntilLengthReachedAnalytical
       (const Vec3& xP, const UnitVec3& tP, const Real& terminatingLength, 
        const GeodesicOptions& options, Geodesic& geod) const override;

    void shootGeodesicInDirectionUntilPlaneHitAnalytical
       (const Vec3& xP, const UnitVec3& tP, const Plane& terminatingPlane, 
        const GeodesicOptions& options, Geodesic& geod) const override;

    void calcGeodesicAnalytical
       (const Vec3& xP, const Vec3& xQ, const Vec3& tPhint, const Vec3& tQhint, 
        Geodesic& geod) const override;

    const Function& getImplicitFunction() const override {
        return function;
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createOBBTree();

    Real                    radius;
    SphereImplicitFunction  function;
};



//==============================================================================
//                                ELLIPSOID IMPL
//==============================================================================
class EllipsoidImplicitFunction : public Function {
public:
    EllipsoidImplicitFunction() : ownerp(0) {}
    EllipsoidImplicitFunction(const ContactGeometry::Ellipsoid::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::Ellipsoid::Impl& owner) {ownerp=&owner;}
    Real calcValue(const Vector& x) const override;
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vector& x) const override;
    int getArgumentSize() const override {return 3;}
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::Ellipsoid::Impl* ownerp;// just a ref.; don't delete
};

class ContactGeometry::Ellipsoid::Impl : public ContactGeometryImpl  {
public:
    explicit Impl(const Vec3& radii)
    :   radii(radii),
        curvatures(Vec3(1/radii[0],1/radii[1],1/radii[2])) 
    {   function.setOwner(*this);
        createOBBTree(); }

    ContactGeometryImpl* clone() const override {return new Impl(radii);}
    const Vec3& getRadii() const {return radii;}
    void setRadii(const Vec3& r) 
    {   radii = r; curvatures = Vec3(1/r[0],1/r[1],1/r[2]); }

    const Vec3& getCurvatures() const {return curvatures;}

    // See below.
    inline Vec3 findPointWithThisUnitNormal(const UnitVec3& n) const;
    inline Vec3 findPointInSameDirection(const Vec3& Q) const;
    inline UnitVec3 findUnitNormalAtPoint(const Vec3& Q) const;

    // Cost is findParaboloidAtPointWithNormal + about 40 flops.
    void findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const
    {   findParaboloidAtPointWithNormal(Q,findUnitNormalAtPoint(Q),X_EP,k); }

    void findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& n,
        Transform& X_EP, Vec2& k) const;

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const override;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    bool isSmooth() const override {return true;}
    bool isConvex() const override {return true;}
    bool isFinite() const override {return true;}

    // The point furthest in this direction is the unique point whose outward
    // normal is this direction.
    Vec3 calcSupportPoint(const UnitVec3& direction) const override {
        return findPointWithThisUnitNormal(direction);
    }
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const override;
    const Function& getImplicitFunction() const override {
        return function;
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createOBBTree();


    Vec3 radii;
    // The curvatures are calculated whenever the radii are set.
    Vec3 curvatures; // (1/radii[0], 1/radii[1], 1/radii[2])
    EllipsoidImplicitFunction function;
};

// Given an ellipsoid and a unit normal direction, find the unique point on the
// ellipsoid whose outward normal matches. The unnormalized normal at a point
// p=[x y z] is the gradient of the ellipsoid's implicit equation there:
// (1)      n(p) = grad(f(p)) = 2*[x/a^2, y/b^2, z/c^2]
// If we had that, we'd have 
// (2)      p=[n[0]*a^2, n[1]*b^2, n[2]*c^2]/2,
// but instead we're given the normalized normal that has been divided
// by the length of n:
// (3)      nn = n/|n| = s * [x/a^2, y/b^2, z/c^2]
// where s = 2/|n|. 
// 
// We can solve for s using the fact that x,y,z must lie on the ellipsoid so
// |x/a,y/b,z/c|=1. Construct the vector
//          v = [nn[0]*a, nn[1]*b, nn[2]*c] = s*[x/a, y/b, z/c]
// Now we have |v|=s. So n/2 = nn/|v| and we can use equation (2) to solve
// for p. Cost is about 40 flops.
inline Vec3 ContactGeometry::Ellipsoid::Impl::
findPointWithThisUnitNormal(const UnitVec3& nn) const {
    const Real& a=radii[0]; const Real& b=radii[1]; const Real& c=radii[2];
    const Vec3 v  = Vec3(nn[0]*a, nn[1]*b, nn[2]*c);
    const Vec3 p  = Vec3( v[0]*a,  v[1]*b,  v[2]*c) / v.norm();
    return p;
}

// Given a point Q=(x,y,z) measured from ellipse center O, find the intersection 
// of the ray d=Q-O with the ellipse surface. This just requires scaling the 
// direction vector d by a factor s so that f(s*d)=0, that is, 
//       s*|x/a y/b z/c|=1  => s = 1/|x/a y/b z/c|
// Cost is about 40 flops.
inline Vec3 ContactGeometry::Ellipsoid::Impl::
findPointInSameDirection(const Vec3& Q) const {
    Real s = 1/Vec3(Q[0]*curvatures[0], 
                    Q[1]*curvatures[1], 
                    Q[2]*curvatures[2]).norm();
    return s*Q;
}

// The implicit equation of the ellipsoid surface is f(x,y,z)=0 where
// f(x,y,z) = (ka x)^2 + (kb y)^2 (kc z)^2 - 1. Points p inside the ellipsoid
// have f(p)<0, outside f(p)>0. f defines a field in space; its positive
// gradient [df/dx df/dy df/dz] points outward. So, given an ellipsoid with 
// principal curvatures ka,kb,kc and a point Q allegedly on the ellipsoid, the
// outward normal (unnormalized) n at that point is
//    n(p) = grad(f(p)) = 2*[ka^2 x, kb^2 y, kc^2 z]
// so the unit norm we're interested in is nn=n/|n| (the "2" drops out).
// If Q is not on the ellipsoid this is equivalent to scaling the ray Q-O
// until it hits the ellipsoid surface at Q'=s*Q, and then reporting the outward
// normal at Q' instead.
// Cost is about 40 flops.
inline UnitVec3 ContactGeometry::Ellipsoid::Impl::
findUnitNormalAtPoint(const Vec3& Q) const {
    const Vec3 kk(square(curvatures[0]), square(curvatures[1]), 
                  square(curvatures[2]));
    return UnitVec3(kk[0]*Q[0], kk[1]*Q[1], kk[2]*Q[2]);
}



//==============================================================================
//                          SMOOTH HEIGHT MAP IMPL
//==============================================================================
class SmoothHeightMapImplicitFunction : public Function {
public:
    SmoothHeightMapImplicitFunction() : ownerp(0) {}
    SmoothHeightMapImplicitFunction
       (const ContactGeometry::SmoothHeightMap::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::SmoothHeightMap::Impl& owner) 
    {   ownerp=&owner; }
    Real calcValue(const Vector& x) const override;
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vector& x) const override;
    int getArgumentSize() const override {return 3;}
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
private:
    // just a reference; don't delete
    const ContactGeometry::SmoothHeightMap::Impl*   ownerp; 
};



class ContactGeometry::SmoothHeightMap::Impl : public ContactGeometryImpl {
public:
    explicit Impl(const BicubicSurface& surface);

    ContactGeometryImpl* clone() const override {
        return new Impl(surface);
    }

    const BicubicSurface& getBicubicSurface() const {return surface;}
    BicubicSurface::PatchHint& updHint() const {return hint;}

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const override;

    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const override;

    void getBoundingSphere(Vec3& center, Real& radius) const override {
        center = boundingSphere.getCenter();
        radius = boundingSphere.getRadius();
    }

    bool isSmooth() const override {return true;}
    bool isConvex() const override {return false;}
    bool isFinite() const override {return true;}

    Vec3 calcSupportPoint(const UnitVec3& direction) const override {
        assert(false);
        return Vec3(NaN);
    }

    // We ignore the z coordinate here and just return the curvature of
    // the unique point at (x,y).
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const override {
        Transform X_SP;
        surface.calcParaboloid(Vec2(point[0],point[1]), hint, X_SP, curvature);
        orientation = X_SP.R();
    }

    const Function& getImplicitFunction() const override 
    {   return implicitFunction; }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createBoundingVolumes();
    // The given OBBNode is assigned this range of patches. If there is
    // more than one patch in the range, it will dole those out to its
    // children recursively until the leaves each have responsibility for
    // one patch. Then we'll deal with that patch, which may have to be
    // subdivided into submission.
    void splitPatches(int x0,int y0, int nx, int ny, 
                      OBBNode& node, int depth,
                      Array_<const Vec3*>* parentControlPoints=0) const;

    // The supplied OBBNode has responsibility for the given subpatch, which
    // may need further subdivision.
    void assignPatch(const Geo::BicubicBezierPatch& patch,
                     OBBNode& node, int depth, 
                     Array_<const Vec3*>* parentControlPoints=0) const;


    BicubicSurface                      surface;
    mutable BicubicSurface::PatchHint   hint;
    Geo::Sphere                         boundingSphere;
    SmoothHeightMapImplicitFunction     implicitFunction;
};





//==============================================================================
//                                BRICK IMPL
//==============================================================================

class ContactGeometry::Brick::Impl : public ContactGeometryImpl {
public:
    explicit Impl(const Vec3& halfLengths) : m_box(halfLengths) {
        createOBBTree(); 
    }

    ContactGeometryImpl* clone() const override {
        return new Impl(getHalfLengths());
    }
    const Vec3& getHalfLengths() const {return m_box.getHalfLengths();}

    void setHalfLengths(const Vec3& halfLengths) {
        m_box.setHalfLengths(halfLengths);
    }

    const Geo::Box& getGeoBox() const {return m_box;}

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const override;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    bool isSmooth() const override {return false;}
    bool isConvex() const override {return true;}
    bool isFinite() const override {return true;}

    Vec3 calcSupportPoint(const UnitVec3& direction) const override {
        return m_box.findSupportPoint(direction);
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createOBBTree();

    Geo::Box        m_box;
};



//==============================================================================
//                            OBB TREE NODE IMPL
//==============================================================================
class OBBTreeNodeImpl {
public:
    OBBTreeNodeImpl() : child1(NULL), child2(NULL), numTriangles(0) {
    }
    OBBTreeNodeImpl(const OBBTreeNodeImpl& copy);
    ~OBBTreeNodeImpl();
    OrientedBoundingBox bounds;
    OBBTreeNodeImpl* child1;
    OBBTreeNodeImpl* child2;
    Array_<int> triangles;
    int numTriangles;
    Vec3 findNearestPoint(const ContactGeometry::TriangleMesh::Impl& mesh, 
                          const Vec3& position, Real cutoff2, Real& distance2, 
                          int& face, Vec2& uv) const;
    bool intersectsRay(const ContactGeometry::TriangleMesh::Impl& mesh, 
                       const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, int& face, Vec2& uv) const;
};



//==============================================================================
//                            TRIANGLE MESH IMPL
//==============================================================================
class ContactGeometry::TriangleMesh::Impl : public ContactGeometryImpl {
public:
    class Edge;
    class Face;
    class Vertex;

    Impl(const ArrayViewConst_<Vec3>& vertexPositions, 
         const ArrayViewConst_<int>& faceIndices, bool smooth);
    Impl(const PolygonalMesh& mesh, bool smooth);
    ContactGeometryImpl* clone() const override {
        return new Impl(*this);
    }

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    Vec3     findPoint(int face, const Vec2& uv) const;
    Vec3     findCentroid(int face) const;
    UnitVec3 findNormalAtPoint(int face, const Vec2& uv) const;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, int& face, 
                          Vec2& uv) const;
    Vec3 findNearestPointToFace(const Vec3& position, int face, Vec2& uv) const;
    void createPolygonalMesh(PolygonalMesh& mesh) const;

    DecorativeGeometry createDecorativeGeometry() const override;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, int& face, Vec2& uv) const;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    bool isSmooth() const override {return false;}
    bool isConvex() const override {return false;}
    bool isFinite() const override {return true;}


    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void init(const Array_<Vec3>& vertexPositions, const Array_<int>& faceIndices);
    void createObbTree(OBBTreeNodeImpl& node, const Array_<int>& faceIndices);
    void splitObbAxis(const Array_<int>& parentIndices, 
                      Array_<int>& child1Indices, 
                      Array_<int>& child2Indices, int axis);
    void findBoundingSphere(Vec3* point[], int p, int b, 
                            Vec3& center, Real& radius);
    friend class ContactGeometry::TriangleMesh;
    friend class OBBTreeNodeImpl;

    Array_<Edge>    edges;
    Array_<Face>    faces;
    Array_<Vertex>  vertices;
    Vec3            boundingSphereCenter;
    Real            boundingSphereRadius;
    OBBTreeNodeImpl obb;
    bool            smooth;
};



//==============================================================================
//                          TriangleMeshImpl EDGE
//==============================================================================
class ContactGeometry::TriangleMesh::Impl::Edge {
public:
    Edge(int vert1, int vert2, int face1, int face2) {
        vertices[0] = vert1;
        vertices[1] = vert2;
        faces[0] = face1;
        faces[1] = face2;
    }
    int     vertices[2];
    int     faces[2];
};



//==============================================================================
//                           TriangleMeshImpl FACE
//==============================================================================
class ContactGeometry::TriangleMesh::Impl::Face {
public:
    Face(int vert1, int vert2, int vert3, 
         const Vec3& normal, Real area) 
    :   normal(normal), area(area) {
        vertices[0] = vert1;
        vertices[1] = vert2;
        vertices[2] = vert3;
        edges[0] = edges[1] = edges[2] = -1; // filled in later
    }
    int         vertices[3];
    int         edges[3];
    UnitVec3    normal;
    Real        area;
};



//==============================================================================
//                          TriangleMeshImpl VERTEX
//==============================================================================
class ContactGeometry::TriangleMesh::Impl::Vertex {
public:
    Vertex(Vec3 pos) : pos(pos), firstEdge(-1) {
    }
    Vec3        pos;
    UnitVec3    normal;
    int         firstEdge;
};



//==============================================================================
//                              TORUS IMPL
//==============================================================================
class TorusImplicitFunction : public Function {
public:
    TorusImplicitFunction() : ownerp(0) {}
    TorusImplicitFunction(const ContactGeometry::Torus::Impl& owner)
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::Torus::Impl& owner) {ownerp=&owner;}
    Real calcValue(const Vector& x) const override;
    Real calcDerivative(const Array_<int>& derivComponents,
                        const Vector& x) const override;
    int getArgumentSize() const override {return 3;}
    int getMaxDerivativeOrder() const override
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::Torus::Impl* ownerp; // just a reference; don't delete
};

class ContactGeometry::Torus::Impl : public ContactGeometryImpl {
public:
    explicit Impl(Real torusRadius, Real tubeRadius) :
    torusRadius(torusRadius), tubeRadius(tubeRadius) {
        function.setOwner(*this);
    }

    ContactGeometryImpl* clone() const override {
        return new Impl(torusRadius, tubeRadius);
    }
    Real getTorusRadius() const {
        return torusRadius;
    }
    void setTorusRadius(Real r) {
        torusRadius = r;
    }
    Real getTubeRadius() const {
        return tubeRadius;
    }
    void setTubeRadius(Real r) {
        tubeRadius = r;
    }

    ContactGeometryTypeId getTypeId() const override {return classTypeId();}

    DecorativeGeometry createDecorativeGeometry() const override;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction,
                       Real& distance, UnitVec3& normal) const override;
    void getBoundingSphere(Vec3& center, Real& radius) const override;

    void createPolygonalMesh(PolygonalMesh& mesh) const;

    bool isSmooth() const override {return true;}
    bool isConvex() const override {return false;}
    bool isFinite() const override {return true;}

    Vec3 calcSupportPoint(const UnitVec3& direction) const override;

    void calcCurvature(const Vec3& point, Vec2& curvature,
                       Rotation& orientation) const override;

    Vec3 findNearestPoint(const Vec3& position, bool& inside,
            UnitVec3& normal) const override;

//    TODO
//    virtual void shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
//            const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const;
//
//    virtual void shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
//            const Plane& terminatingPlane, const GeodesicOptions& options,
//            Geodesic& geod) const;
//
//    virtual void calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
//                const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const;

    const Function& getImplicitFunction() const override {
        return function;
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id =
            createNewContactGeometryTypeId();
        return id;
    }
private:
    Real                    torusRadius;
    Real                    tubeRadius;
    TorusImplicitFunction   function;
};


} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_GEOMETRY_IMPL_H_
