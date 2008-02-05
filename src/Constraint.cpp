/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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

/**@file
 *
 * Private implementation of MobilizedBody, and its included subclasses which
 * represent the built-in mobilizer types.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintRep.h"
#include "SimbodyMatterSubsystemRep.h"

#ifdef USE_OLD_CONSTRAINTS
    #include "LengthConstraints.h"
#endif

namespace SimTK {


    ////////////////
    // CONSTRAINT //
    ////////////////

bool Constraint::isEmptyHandle() const {return rep==0;}
bool Constraint::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}

void Constraint::disown(Constraint& newOwnerHandle) {
    SimTK_ASSERT_ALWAYS(rep && rep->myHandle==this,
        "disown() not allowed for an empty or non-owner Constraint handle.");
    SimTK_ASSERT_ALWAYS(!newOwnerHandle.rep,
        "disown() can only transfer ownership to an empty Constraint handle.");

    newOwnerHandle.setRep(*rep);
    rep->setMyHandle(newOwnerHandle);
}

Constraint::~Constraint() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

// Make this Constraint a non-owner handle referring to the same
// object as the source.
Constraint::Constraint(Constraint& src) : rep(src.rep) {
}

// Make this empty or non-owner handle refer to the same object
// as the source. This is illegal if the current handle is an
// owner.
Constraint& Constraint::operator=(Constraint& src) {
    if (&src != this) {
        SimTK_ASSERT_ALWAYS(!(rep && rep->myHandle==this),
            "You can't reassign the owner handle of a Constraint.");
        rep = src.rep;
    }
    return *this;
}

const SimbodyMatterSubsystem& Constraint::getMatterSubsystem() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getMatterSubsystem() called on a Constraint that is not part of a subsystem.");
    return getRep().getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}

ConstraintIndex Constraint::getConstraintIndex() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getConstraintIndex() called on a Constraint that is not part of a subsystem.");
    return rep->myConstraintIndex;
}

SimbodyMatterSubsystem& Constraint::updMatterSubsystem() {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "updMatterSubsystem() called on a Constraint that is not part of a subsystem.");
    return updRep().updMyMatterSubsystemRep().updMySimbodyMatterSubsystemHandle();
}

bool Constraint::isInSubsystem() const {
    return getRep().isInSubsystem();
}

bool Constraint::isInSameSubsystem(const MobilizedBody& body) const {
    return getRep().isInSameSubsystem(body);
}

int Constraint::getNumConstrainedBodies() const {
    return getRep().getNumConstrainedBodies();
}

int Constraint::getNumConstrainedMobilities(const State& s) const {
    return getRep().getNumConstrainedMobilities(s);
}

int Constraint::getNumConstrainedMobilities(const State& s, ConstrainedBodyIndex B) const {
    return getRep().getNumConstrainedMobilities(s,B);
}

int Constraint::getConstrainedMobilityIndex(const State& s, ConstrainedBodyIndex B, int which) const {
    return getRep().getConstrainedMobilityIndex(s,B,which);
}

const MobilizedBody& Constraint::getConstrainedMobilizedBody(ConstrainedBodyIndex B) const {
    return getRep().getConstrainedMobilizedBody(B);
}
const MobilizedBody& Constraint::getAncestorMobilizedBody() const {
    return getRep().getAncestorMobilizedBody();
}

const SimbodyMatterSubsystem::Subtree& Constraint::getSubtree() const {
    assert(getRep().subsystemTopologyHasBeenRealized());
    return getRep().mySubtree;
}

// Find out how many holonomic (position), nonholonomic (velocity),
// and acceleration-only constraint equations are generated by this Constraint.
void Constraint::getNumConstraintEquations(const State& s, int& mp, int& mv, int& ma) const {
	getRep().getNumConstraintEquations(s,mp,mv,ma);
}

Vector Constraint::getPositionError(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	Vector perr(mp);
	if (mp) getRep().getPositionErrors(s, mp, &perr[0]);
	return perr;
}

Vector Constraint::getVelocityError(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	Vector pverr(mp+mv);
	if (mp+mv) getRep().getVelocityErrors(s, mp+mv, &pverr[0]);
	return pverr;
}

Vector Constraint::getAccelerationError(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	Vector pvaerr(mp+mv+ma);
	if (mp+mv+ma) getRep().getAccelerationErrors(s, mp+mv+ma, &pvaerr[0]);
	return pvaerr;
}

Matrix Constraint::calcPositionConstraintMatrixP(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);

	Matrix P(mp, nu);
	if (mp && nu) {
		Vector  pverr0(mp), pverr(mp); // we're interested in the first mp of these
		State   tmp = s;      // don't change s

		matter.updU(tmp) = 0;	// first calculate the bias term -c(t,q)
		system.realize(tmp, Stage::Velocity);
		pverr0 = getVelocityError(tmp)(0,mp);

		// Now calculate sensitivity of d(perr)/dt=Pu-c(t,q) to each u in turn.
		for (int j=0; j<nu; ++j) {
			matter.updU(tmp)[j] = 1;
		    system.realize(tmp, Stage::Velocity);
			pverr = getVelocityError(tmp)(0,mp);
			matter.updU(tmp)[j] = 0;
			P(j) = pverr - pverr0;
		}
	}
	return P;
}

Matrix Constraint::calcPositionConstraintMatrixPt(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nb = matter.getNBodies();

	const int ncb = getNumConstrainedBodies();

	Matrix Pt(nu, mp);
	if (mp && nu) {
		const ConstraintRep& rep = getRep();
		Vector_<SpatialVec> bodyForcesInA(ncb);
		Vector              mobilityForces(nu); // TODO should be n participating u's

		Vector_<SpatialVec> bodyForcesInG(nb);
		bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

		// Calculate Pt*lambda with each lambda set to 1 in turn.
		Vector lambda(mp);
		lambda = 0;
		for (int i=0; i<mp; ++i) {
			lambda[i] = 1;
			bodyForcesInA = SpatialVec(Vec3(0), Vec3(0));
			mobilityForces = 0;
			rep.applyPositionConstraintForces(s, mp, &lambda[0], bodyForcesInA, mobilityForces);
			for (ConstrainedBodyIndex cb(0); cb < ncb; ++cb) {
				bodyForcesInG[rep.getMobilizedBodyIndexOfConstrainedBody(cb)] =
					R_GA*bodyForcesInA[cb];
			}
			lambda[i] = 0;

			rep.getMyMatterSubsystem().calcInternalGradientFromSpatial(s,bodyForcesInG,Pt(i));
            //TODO: must unpack and add in mobilityForces
		}

	}
	return Pt;
}

// Calculate the constraint matrix V= partial(verr)/partial(u) for just
// the nonholonomic constraints. For this we use the acceleration error
// equations obtained from differentiation of the nonholonomic constraints:
//    verr_dot = V udot - b(t,q,u)
//
Matrix Constraint::calcVelocityConstraintMatrixV(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s); // same as nudot

	Matrix V(mv, nu);
	if (mv && nu) {
		Vector  vaerr0(mv), vaerr(mv); // we're interested in the middle mv of these (mp-mv-ma)

		Vector udot(nu);
		udot = 0;

		// Calculate the bias term -b(t,q,u)
		vaerr0 = calcAccelerationErrorFromUDot(s, udot)(mp, mv);

		// Now calculate sensitivity of d(verr)/dt=Vudot-b(t,q,u) to each udot in turn.
		for (int j=0; j<nu; ++j) {
			udot[j] = 1;
		    vaerr = calcAccelerationErrorFromUDot(s, udot)(mp, mv);
			udot[j] = 0;
			V(j) = vaerr - vaerr0;
		}
	}
	return V;
}

Matrix Constraint::calcVelocityConstraintMatrixVt(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nb = matter.getNBodies();

	const int ncb = getNumConstrainedBodies();

	Matrix Vt(nu, mv);
	if (mv && nu) {
		const ConstraintRep& rep = getRep();
		Vector_<SpatialVec> bodyForcesInA(ncb);
		Vector              mobilityForces(nu); // TODO should be n participating u's

		Vector_<SpatialVec> bodyForcesInG(nb);

		bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

		// Calculate Vt*lambda with each lambda set to 1 in turn.
		Vector lambda(mv);
		lambda = 0;
		for (int i=0; i<mv; ++i) {
			lambda[i] = 1;
			bodyForcesInA = SpatialVec(Vec3(0), Vec3(0));
			mobilityForces = 0;
			rep.applyVelocityConstraintForces(s, mv, &lambda[0], bodyForcesInA, mobilityForces);
			for (ConstrainedBodyIndex cb(0); cb < ncb; ++cb) {
				bodyForcesInG[rep.getMobilizedBodyIndexOfConstrainedBody(cb)] =
					R_GA*bodyForcesInA[cb];
			}
			lambda[i] = 0;

			rep.getMyMatterSubsystem().calcInternalGradientFromSpatial(s,bodyForcesInG,Vt(i));
            //TODO: must unpack and add in mobilityForces
		}

	}
	return Vt;
}

// Calculate the constraint matrix A= partial(aerr)/partial(udot) for just
// the acceleration-only constraints. For this we use the acceleration error
// equations directly because we're *requiring* that acceleration-only constraints
// are linear in the udots.
//    aerr = A udot - b(t,q,u)
//
Matrix Constraint::calcAccelerationConstraintMatrixA(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nudot = matter.getNU(s); // same as nudot

	Matrix A(ma, nudot);
	if (ma && nudot) {
		Vector  aerr0(ma), aerr(ma); // we're interested in the final ma of these (mp-mv-ma)

		Vector udot(nudot);
		udot = 0;

		// Calculate the bias term -b(t,q,u)
		aerr0 = calcAccelerationErrorFromUDot(s, udot)(mp+mv, ma);

		// Now calculate sensitivity of d(aerr)/dt=Audot-b(t,q,u) to each udot in turn.
		for (int j=0; j<nudot; ++j) {
			udot[j] = 1;
		    aerr = calcAccelerationErrorFromUDot(s, udot)(mp+mv, ma);
			udot[j] = 0;
			A(j) = aerr - aerr0;
		}
	}
	return A;
}

Matrix Constraint::calcAccelerationConstraintMatrixAt(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nb = matter.getNBodies();

	const int ncb = getNumConstrainedBodies();

	Matrix At(nu, ma);
	if (ma && nu) {
		const ConstraintRep& rep = getRep();
		Vector_<SpatialVec> bodyForcesInA(ncb);
		Vector              mobilityForces(nu); // TODO should be n participating u's

		Vector_<SpatialVec> bodyForcesInG(nb);
		bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

		// Calculate At*lambda with each lambda set to 1 in turn.
		Vector lambda(ma);
		lambda = 0;
		for (int i=0; i<ma; ++i) {
			lambda[i] = 1;
			bodyForcesInA = SpatialVec(Vec3(0), Vec3(0));
			mobilityForces = 0;
			rep.applyAccelerationConstraintForces(s, ma, &lambda[0], bodyForcesInA, mobilityForces);
			for (ConstrainedBodyIndex cb(0); cb < ncb; ++cb) {
				bodyForcesInG[rep.getMobilizedBodyIndexOfConstrainedBody(cb)] =
					R_GA*bodyForcesInA[cb];
			}
			lambda[i] = 0;

			rep.getMyMatterSubsystem().calcInternalGradientFromSpatial(s,bodyForcesInG,At(i));
            //TODO: must unpack and add in mobilityForces
		}

	}
	return At;
}

Matrix Constraint::calcPositionConstraintMatrixPQInverse(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquations(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nq = matter.getNQ(s);

	const Matrix P = calcPositionConstraintMatrixP(s);
	assert(P.nrow()==mp && P.ncol()==nu);

	Matrix PQInv(mp, nq); // = P*Q^-1
	if (mp && nq) {
		// The routine below calculates qlikeRow = ulikeRow * Q^-1 which is what
		// we need but it actually works with Vectors (transpose of RowVectors)
		// and at the moment they have to be contiguous so we'll have to copy.
		Vector uin(nu);
		Vector qout(nq);
		for (int i=0; i < mp; ++i) {
			uin = ~P[i];
			matter.multiplyByQMatrixInverse(s, true, uin, qout);
			PQInv[i] = ~qout;
		}
	}
	return PQInv;
}

void Constraint::calcConstraintForcesFromMultipliers(
    const State&         s,
    const Vector&        lambda,
    Vector_<SpatialVec>& bodyForcesInA,
    Vector&              mobilityForces) const
{
    int mp, mv, ma;
    getNumConstraintEquations(s, mp, mv, ma);
    assert(lambda.size() == mp+mv+ma);
    assert(lambda.hasContiguousData());

    getRep().calcConstraintForcesFromMultipliers(s,mp,mv,ma,&lambda[0],bodyForcesInA,mobilityForces);
}


    /////////////////////
    // CONSTRAINT::ROD //
    /////////////////////

Constraint::Rod::Rod(MobilizedBody& body1, MobilizedBody& body2, Real defaultRodLength)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    rep = new RodRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B1 = updRep().addConstrainedBody(body1);
    updRep().B2 = updRep().addConstrainedBody(body2);

    updRep().defaultRodLength = defaultRodLength;
}

Constraint::Rod::Rod(MobilizedBody& body1, const Vec3& point1,
                     MobilizedBody& body2, const Vec3& point2, Real defaultRodLength)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    rep = new RodRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B1 = updRep().addConstrainedBody(body1);
    updRep().B2 = updRep().addConstrainedBody(body2);

    updRep().defaultPoint1 = point1;
    updRep().defaultPoint2 = point2;
    updRep().defaultRodLength = defaultRodLength;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody1(const Vec3& p1) {
    updRep().defaultPoint1 = p1;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody2(const Vec3& p2) {
    updRep().defaultPoint2 = p2;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultRodLength(Real length) {
    updRep().defaultRodLength = length;
    return *this;
}


MobilizedBodyIndex Constraint::Rod::getBody1MobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B1);
}
MobilizedBodyIndex Constraint::Rod::getBody2MobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B1);
}
const Vec3& Constraint::Rod::getDefaultPointOnBody1() const {
    return getRep().defaultPoint1;
}
const Vec3& Constraint::Rod::getDefaultPointOnBody2() const {
    return getRep().defaultPoint2;
}
Real Constraint::Rod::getDefaultRodLength() const {
    return getRep().defaultRodLength;
}



    // Rod bookkeeping //

bool Constraint::Rod::isInstanceOf(const Constraint& s) {
    return RodRep::isA(s.getRep());
}
const Constraint::Rod& Constraint::Rod::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Rod&>(s);
}
Constraint::Rod& Constraint::Rod::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Rod&>(s);
}
const Constraint::Rod::RodRep& Constraint::Rod::getRep() const {
    return dynamic_cast<const RodRep&>(*rep);
}
Constraint::Rod::RodRep& Constraint::Rod::updRep() {
    return dynamic_cast<RodRep&>(*rep);
}
    // RodRep

void Constraint::Rod::RodRep::realizeTopologyVirtual(State& s) const { 
#ifdef USE_OLD_CONSTRAINTS
    SimbodyMatterSubsystemRep& matter = 
        const_cast<RodRep*>(this)->updMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B1);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(B2);
    const RigidBodyNode& rbn1 = matter.getRigidBodyNode(mobilizedBody1);
    const RigidBodyNode& rbn2 = matter.getRigidBodyNode(mobilizedBody2);
    const RBStation s1(rbn1, defaultPoint1);
    const RBStation s2(rbn2, defaultPoint2);
    matter.addOneDistanceConstraintEquation(s1,s2,defaultRodLength);
#endif
}

void Constraint::Rod::RodRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the endpoint artwork until we know the end point stations,
    // which could be as late as Stage::Instance.
    if (stage == Stage::Instance && pointRadius != 0 && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        // TODO: point stations and rod length should be instance-stage data 
        // from State rather than topological data
        const MobilizedBodyIndex body1 = getMobilizedBodyIndexOfConstrainedBody(B1);
        const MobilizedBodyIndex body2 = getMobilizedBodyIndexOfConstrainedBody(B2);

        const Real useRadius = pointRadius > 0 ? pointRadius 
            : std::max(defaultRodLength, .1) * 0.02; // 2% of the length by default

        // Draw a blue mesh sphere at the first point.
        geom.push_back(DecorativeSphere(useRadius)
                            .setColor(Blue)
                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                            .setResolution(0.5)
                            .setBodyId(body1)
                            .setTransform(defaultPoint1));

        // On the follower body draw an purple mesh sphere at the point radius.
        geom.push_back(DecorativeSphere(useRadius)
                            .setColor(Purple)
                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                            .setResolution(0.5)
                            .setBodyId(body2)
                            .setTransform(defaultPoint2));
    }

    // We can't generate the line artwork until we know the two end point locations,
    // which isn't until Position stage since the ends are on different bodies.
    if (stage == Stage::Position) {
        // TODO: point stations and rod length should be instance-stage data 
        // from State rather than topological data

        const Vec3 p_GP1 = getConstrainedMobilizedBody(B1)
                              .locateBodyPointOnGround(s, defaultPoint1);
        const Vec3 p_GP2 = getConstrainedMobilizedBody(B2)
                              .locateBodyPointOnGround(s, defaultPoint2);

        const Vec3 p_P1P2 = p_GP2 - p_GP1;
        const Real d = p_P1P2.norm();

        if (d >= SignificantReal) {
            const Vec3 endPoint = p_GP1 + defaultRodLength * p_P1P2/d;
            geom.push_back(DecorativeLine(p_GP1, endPoint)
                                            .setColor(Black)
                                            .setLineThickness(3)
                                            .setBodyId(GroundIndex));
        }
    }

}

    ////////////////////////////////
    // CONSTRAINT::POINT IN PLANE //
    ////////////////////////////////

Constraint::PointInPlane::PointInPlane
   (MobilizedBody& planeBody,    const UnitVec3& defPlaneNormal, Real defPlaneHeight,
    MobilizedBody& followerBody, const Vec3&     defFollowerPoint)
{
    SimTK_ASSERT_ALWAYS(planeBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::PointInPlane(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(planeBody.isInSameSubsystem(followerBody),
        "Constraint::PointInPlane(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    rep = new PointInPlaneRep(); rep->setMyHandle(*this);
    planeBody.updMatterSubsystem().adoptConstraint(*this);

    updRep().planeBody    = updRep().addConstrainedBody(planeBody);
    updRep().followerBody = updRep().addConstrainedBody(followerBody);
    updRep().defaultPlaneNormal   = defPlaneNormal;
    updRep().defaultPlaneHeight   = defPlaneHeight;
    updRep().defaultFollowerPoint = defFollowerPoint;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultPlaneNormal(const UnitVec3& n) {
    getRep().invalidateTopologyCache();
    updRep().defaultPlaneNormal = n;
    return *this;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultPlaneHeight(Real h) {
    getRep().invalidateTopologyCache();
    updRep().defaultPlaneHeight = h;
    return *this;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultFollowerPoint(const Vec3& p) {
    getRep().invalidateTopologyCache();
    updRep().defaultFollowerPoint = p;
    return *this;
}

MobilizedBodyIndex Constraint::PointInPlane::getPlaneMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().planeBody);
}
MobilizedBodyIndex Constraint::PointInPlane::getFollowerMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().followerBody);
}
const UnitVec3& Constraint::PointInPlane::getDefaultPlaneNormal() const {
    return getRep().defaultPlaneNormal;
}
Real Constraint::PointInPlane::getDefaultPlaneHeight() const {
    return getRep().defaultPlaneHeight;
}
const Vec3& Constraint::PointInPlane::getDefaultFollowerPoint() const {
    return getRep().defaultFollowerPoint;
}

Constraint::PointInPlane& Constraint::PointInPlane::setPlaneDisplayHalfWidth(Real h) {
    updRep().setPlaneDisplayHalfWidth(h);
    return *this;
}
Constraint::PointInPlane& Constraint::PointInPlane::setPointDisplayRadius(Real r) {
    updRep().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::PointInPlane::getPlaneDisplayHalfWidth() const {
    return getRep().getPlaneDisplayHalfWidth();
}

Real Constraint::PointInPlane::getPointDisplayRadius() const {
    return getRep().getPointDisplayRadius();
}

    // PointInPlane bookkeeping //

bool Constraint::PointInPlane::isInstanceOf(const Constraint& s) {
    return PointInPlaneRep::isA(s.getRep());
}
const Constraint::PointInPlane& Constraint::PointInPlane::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const PointInPlane&>(s);
}
Constraint::PointInPlane& Constraint::PointInPlane::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<PointInPlane&>(s);
}
const Constraint::PointInPlane::PointInPlaneRep& Constraint::PointInPlane::getRep() const {
    return dynamic_cast<const PointInPlaneRep&>(*rep);
}

Constraint::PointInPlane::PointInPlaneRep& Constraint::PointInPlane::updRep() {
    return dynamic_cast<PointInPlaneRep&>(*rep);
}

    // PointInPlaneRep

void Constraint::PointInPlane::PointInPlaneRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the normal, height, and follower
    // point location, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        // TODO: should be instance-stage data from State rather than topological data
        // This makes z axis point along plane normal
        const Transform X_B1(Rotation(defaultPlaneNormal,ZAxis), defaultPlaneHeight*defaultPlaneNormal);
        const Transform X_B2(Rotation(), defaultFollowerPoint);

        const MobilizedBodyIndex planeMBId = getMobilizedBodyIndexOfConstrainedBody(planeBody);
        const MobilizedBodyIndex followerMBId = getMobilizedBodyIndexOfConstrainedBody(followerBody);

        if (planeHalfWidth > 0 && pointRadius > 0) {
            // On the inboard body, draw a gray transparent rectangle, outlined in black lines.
            geom.push_back(DecorativeBrick(Vec3(planeHalfWidth,planeHalfWidth,pointRadius/2))
                                                .setColor(Gray)
                                                .setRepresentation(DecorativeGeometry::DrawSurface)
                                                .setOpacity(0.3)
                                                .setBodyId(planeMBId)
                                                .setTransform(X_B1));
            geom.push_back(DecorativeBrick(Vec3(planeHalfWidth,planeHalfWidth,pointRadius/2))
                                                .setColor(Black)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setBodyId(planeMBId)
                                                .setTransform(X_B1));

            // On the follower body draw an orange mesh sphere at the point radius.
            geom.push_back(DecorativeSphere(pointRadius)
                                                .setColor(Orange)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setResolution(0.5)
                                                .setBodyId(followerMBId)
                                                .setTransform(X_B2));
        }
    }
}



    ///////////////////////////////
    // CONSTRAINT::POINT ON LINE //
    ///////////////////////////////

Constraint::PointOnLine::PointOnLine
   (MobilizedBody& lineBody,     const UnitVec3& defLineDirection, const Vec3& defPointOnLine,
    MobilizedBody& followerBody, const Vec3&     defFollowerPoint)
{
    SimTK_ASSERT_ALWAYS(lineBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::PointOnLine(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(lineBody.isInSameSubsystem(followerBody),
        "Constraint::PointOnLine(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    rep = new PointOnLineRep(); rep->setMyHandle(*this);
    lineBody.updMatterSubsystem().adoptConstraint(*this);

    updRep().lineBody     = updRep().addConstrainedBody(lineBody);
    updRep().followerBody = updRep().addConstrainedBody(followerBody);
    updRep().defaultLineDirection = defLineDirection;
    updRep().defaultPointOnLine   = defPointOnLine;
    updRep().defaultFollowerPoint = defFollowerPoint;
}

Constraint::PointOnLine& Constraint::PointOnLine::setDefaultLineDirection(const UnitVec3& z) {
    getRep().invalidateTopologyCache();
    updRep().defaultLineDirection = z;
    return *this;
}

Constraint::PointOnLine& Constraint::PointOnLine::setDefaultPointOnLine(const Vec3& P) {
    getRep().invalidateTopologyCache();
    updRep().defaultPointOnLine = P;
    return *this;
}

Constraint::PointOnLine& Constraint::PointOnLine::setDefaultFollowerPoint(const Vec3& S) {
    getRep().invalidateTopologyCache();
    updRep().defaultFollowerPoint = S;
    return *this;
}

MobilizedBodyIndex Constraint::PointOnLine::getLineMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().lineBody);
}
MobilizedBodyIndex Constraint::PointOnLine::getFollowerMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().followerBody);
}
const UnitVec3& Constraint::PointOnLine::getDefaultLineDirection() const {
    return getRep().defaultLineDirection;
}
const Vec3& Constraint::PointOnLine::getDefaultPointOnLine() const {
    return getRep().defaultPointOnLine;
}
const Vec3& Constraint::PointOnLine::getDefaultFollowerPoint() const {
    return getRep().defaultFollowerPoint;
}

Constraint::PointOnLine& Constraint::PointOnLine::setLineDisplayHalfLength(Real h) {
    updRep().setLineDisplayHalfLength(h);
    return *this;
}
Constraint::PointOnLine& Constraint::PointOnLine::setPointDisplayRadius(Real r) {
    updRep().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::PointOnLine::getLineDisplayHalfLength() const {
    return getRep().getLineDisplayHalfLength();
}

Real Constraint::PointOnLine::getPointDisplayRadius() const {
    return getRep().getPointDisplayRadius();
}

    // PointOnLine bookkeeping //

bool Constraint::PointOnLine::isInstanceOf(const Constraint& s) {
    return PointOnLineRep::isA(s.getRep());
}
const Constraint::PointOnLine& Constraint::PointOnLine::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const PointOnLine&>(s);
}
Constraint::PointOnLine& Constraint::PointOnLine::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<PointOnLine&>(s);
}
const Constraint::PointOnLine::PointOnLineRep& Constraint::PointOnLine::getRep() const {
    return dynamic_cast<const PointOnLineRep&>(*rep);
}

Constraint::PointOnLine::PointOnLineRep& Constraint::PointOnLine::updRep() {
    return dynamic_cast<PointOnLineRep&>(*rep);
}

    // PointOnLineRep

void Constraint::PointOnLine::PointOnLineRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the direction, point on line, and follower
    // point location, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        // TODO: should be instance-stage data from State rather than topological data
        // This makes z axis point along line
        const Transform X_B1(Rotation(defaultLineDirection,ZAxis), defaultPointOnLine);
        const Transform X_B2(Rotation(), defaultFollowerPoint);

        const MobilizedBodyIndex lineMBId = getMobilizedBodyIndexOfConstrainedBody(lineBody);
        const MobilizedBodyIndex followerMBId = getMobilizedBodyIndexOfConstrainedBody(followerBody);

        if (lineHalfLength > 0 && pointRadius > 0) {
            // On the line body, draw a black line centered at the point-on-line.
            geom.push_back(DecorativeLine(Vec3(0,0,-lineHalfLength), Vec3(0,0,lineHalfLength))
                                                .setColor(Black)
                                                .setLineThickness(3)
                                                .setBodyId(lineMBId)
                                                .setTransform(X_B1));

            // On the line body draw a blue mesh sphere at the line center.
            geom.push_back(DecorativeSphere(pointRadius)
                                                .setColor(Blue)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setResolution(0.5)
                                                .setBodyId(lineMBId)
                                                .setTransform(X_B1));

            // On the follower body draw an orange mesh sphere at the point radius.
            geom.push_back(DecorativeSphere(pointRadius)
                                                .setColor(Orange)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setResolution(0.5)
                                                .setBodyId(followerMBId)
                                                .setTransform(X_B2));
        }
    }
}



    ////////////////////////////////
    // CONSTRAINT::CONSTANT ANGLE //
    ////////////////////////////////

Constraint::ConstantAngle::ConstantAngle
   (MobilizedBody& baseBody,     const UnitVec3& defaultAxisOnB,
    MobilizedBody& followerBody, const UnitVec3& defaultAxisOnF,
    Real angle)
{
    SimTK_ASSERT_ALWAYS(baseBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::ConstantAngle(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(baseBody.isInSameSubsystem(followerBody),
        "Constraint::ConstantAngle(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    rep = new ConstantAngleRep(); rep->setMyHandle(*this);
    baseBody.updMatterSubsystem().adoptConstraint(*this);

    updRep().B = updRep().addConstrainedBody(baseBody);
    updRep().F = updRep().addConstrainedBody(followerBody);
    updRep().defaultAxisB = defaultAxisOnB;
    updRep().defaultAxisF = defaultAxisOnF;
    updRep().defaultAngle = angle;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setDefaultBaseAxis(const UnitVec3& a) {
    getRep().invalidateTopologyCache();
    updRep().defaultAxisB = a;
    return *this;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setDefaultFollowerAxis(const UnitVec3& a) {
    getRep().invalidateTopologyCache();
    updRep().defaultAxisF = a;
    return *this;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setDefaultAngle(Real t) {
    getRep().invalidateTopologyCache();
    updRep().defaultAngle = t;
    return *this;
}

MobilizedBodyIndex Constraint::ConstantAngle::getBaseMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B);
}
MobilizedBodyIndex Constraint::ConstantAngle::getFollowerMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().F);
}
const UnitVec3& Constraint::ConstantAngle::getDefaultBaseAxis() const {
    return getRep().defaultAxisB;
}
const UnitVec3& Constraint::ConstantAngle::getDefaultFollowerAxis() const {
    return getRep().defaultAxisF;
}
Real Constraint::ConstantAngle::getDefaultAngle() const {
    return getRep().defaultAngle;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setAxisDisplayLength(Real l) {
    updRep().axisLength = l;
    return *this;
}
Constraint::ConstantAngle& Constraint::ConstantAngle::setAxisDisplayWidth(Real w) {
    updRep().axisThickness = w;
    return *this;
}

Real Constraint::ConstantAngle::getAxisDisplayLength() const {
    return getRep().axisLength;
}

Real Constraint::ConstantAngle::getAxisDisplayWidth() const {
    return getRep().axisThickness;
}

    // ConstantAngle bookkeeping //

bool Constraint::ConstantAngle::isInstanceOf(const Constraint& s) {
    return ConstantAngleRep::isA(s.getRep());
}
const Constraint::ConstantAngle& Constraint::ConstantAngle::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const ConstantAngle&>(s);
}
Constraint::ConstantAngle& Constraint::ConstantAngle::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<ConstantAngle&>(s);
}
const Constraint::ConstantAngle::ConstantAngleRep& Constraint::ConstantAngle::getRep() const {
    return dynamic_cast<const ConstantAngleRep&>(*rep);
}

Constraint::ConstantAngle::ConstantAngleRep& Constraint::ConstantAngle::updRep() {
    return dynamic_cast<ConstantAngleRep&>(*rep);
}

    // ConstantAngleRep

void Constraint::ConstantAngle::ConstantAngleRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the normal, height, and follower
    // point location, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        //TODO
    }
}

    //////////////////////
    // CONSTRAINT::BALL //
    //////////////////////

Constraint::Ball::Ball(MobilizedBody& body1, MobilizedBody& body2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Ball(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Ball(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new BallRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B1 = updRep().addConstrainedBody(body1);
    updRep().B2 = updRep().addConstrainedBody(body2);
}

Constraint::Ball::Ball(MobilizedBody& body1, const Vec3& point1,
                       MobilizedBody& body2, const Vec3& point2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Ball(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Ball(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new BallRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B1 = updRep().addConstrainedBody(body1);
    updRep().B2 = updRep().addConstrainedBody(body2);

    updRep().defaultPoint1 = point1;
    updRep().defaultPoint2 = point2;
}

Constraint::Ball& Constraint::Ball::setDefaultPointOnBody1(const Vec3& p1) {
    getRep().invalidateTopologyCache();
    updRep().defaultPoint1 = p1;
    return *this;
}

Constraint::Ball& Constraint::Ball::setDefaultPointOnBody2(const Vec3& p2) {
    getRep().invalidateTopologyCache();
    updRep().defaultPoint2 = p2;
    return *this;
}

MobilizedBodyIndex Constraint::Ball::getBody1MobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B1);
}
MobilizedBodyIndex Constraint::Ball::getBody2MobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B2);
}
const Vec3& Constraint::Ball::getDefaultPointOnBody1() const {
    return getRep().defaultPoint1;
}
const Vec3& Constraint::Ball::getDefaultPointOnBody2() const {
    return getRep().defaultPoint2;
}

Constraint::Ball& Constraint::Ball::setDefaultRadius(Real r) {
    getRep().invalidateTopologyCache();
    updRep().setDefaultRadius(r);
    return *this;
}

Real Constraint::Ball::getDefaultRadius() const {
    return getRep().getDefaultRadius();
}


    // Ball bookkeeping //

bool Constraint::Ball::isInstanceOf(const Constraint& s) {
    return BallRep::isA(s.getRep());
}
const Constraint::Ball& Constraint::Ball::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Ball&>(s);
}
Constraint::Ball& Constraint::Ball::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Ball&>(s);
}
const Constraint::Ball::BallRep& Constraint::Ball::getRep() const {
    return dynamic_cast<const BallRep&>(*rep);
}

Constraint::Ball::BallRep& Constraint::Ball::updRep() {
    return dynamic_cast<BallRep&>(*rep);
}

    // BallRep

void Constraint::Ball::BallRep::realizeTopologyVirtual(State& s) const { 
#ifdef USE_OLD_CONSTRAINTS
    SimbodyMatterSubsystemRep& matter = 
        const_cast<BallRep*>(this)->updMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B1);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(B2);
    const RigidBodyNode& rbn1 = matter.getRigidBodyNode(mobilizedBody1);
    const RigidBodyNode& rbn2 = matter.getRigidBodyNode(mobilizedBody2);
    const RBStation s1x(rbn1, defaultPoint1+Vec3(1,0,0));
    const RBStation s1y(rbn1, defaultPoint1+Vec3(0,1,0));       
    const RBStation s1z(rbn1, defaultPoint1+Vec3(0,0,1));
    const RBStation s2(rbn2, defaultPoint2);
    matter.addOneDistanceConstraintEquation(s1x,s2,1.);
    matter.addOneDistanceConstraintEquation(s1y,s2,1.);
    matter.addOneDistanceConstraintEquation(s1z,s2,1.);
#endif
}

void Constraint::Ball::BallRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the ball until we know the radius, and we can't place
    // the geometry on the body until we know the body1 and body2 point
    // placements on the bodies, which might not be until Instance stage.
    if (stage == Stage::Instance && getDefaultRadius() > 0 && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Transform X_B1(Rotation(), defaultPoint1); // should be point from State
        const Transform X_B2(Rotation(), defaultPoint2);

        // On the inboard body, draw a solid sphere and a wireframe one attached to it for
        // easier visualization of its rotation. These are at about 90% of the radius.
        geom.push_back(DecorativeSphere(0.92*getDefaultRadius())
                                            .setColor(Gray)
                                            .setRepresentation(DecorativeGeometry::DrawSurface)
                                            .setOpacity(0.5)
                                            .setResolution(0.75)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B1))
                                            .setTransform(X_B1));
        geom.push_back(DecorativeSphere(0.90*getDefaultRadius())
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B1))
            .setTransform(X_B1));

        // On the outboard body draw an orange mesh sphere at the ball radius.
        geom.push_back(DecorativeSphere(getDefaultRadius())
                                            .setColor(Orange)
                                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                                            .setOpacity(0.5)
                                            .setResolution(0.5)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B2))
                                            .setTransform(X_B2));
    }
}

    //////////////////////////////////////
    // CONSTRAINT::CONSTANT ORIENTATION //
    //////////////////////////////////////

Constraint::ConstantOrientation::ConstantOrientation
   (MobilizedBody& baseBody,     const Rotation& defaultFrameOnB,
    MobilizedBody& followerBody, const Rotation& defaultFrameOnF)
{
    SimTK_ASSERT_ALWAYS(baseBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::ConstantOrientation(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(baseBody.isInSameSubsystem(followerBody),
        "Constraint::ConstantOrientation(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    rep = new ConstantOrientationRep(); rep->setMyHandle(*this);
    baseBody.updMatterSubsystem().adoptConstraint(*this);

    updRep().B = updRep().addConstrainedBody(baseBody);
    updRep().F = updRep().addConstrainedBody(followerBody);
    updRep().defaultRB = defaultFrameOnB;
    updRep().defaultRF = defaultFrameOnF;
}

Constraint::ConstantOrientation& Constraint::ConstantOrientation::setDefaultBaseRotation(const Rotation& R) {
    getRep().invalidateTopologyCache();
    updRep().defaultRB = R;
    return *this;
}

Constraint::ConstantOrientation& Constraint::ConstantOrientation::setDefaultFollowerRotation(const Rotation& R) {
    getRep().invalidateTopologyCache();
    updRep().defaultRF = R;
    return *this;
}


MobilizedBodyIndex Constraint::ConstantOrientation::getBaseMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B);
}
MobilizedBodyIndex Constraint::ConstantOrientation::getFollowerMobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().F);
}
const Rotation& Constraint::ConstantOrientation::getDefaultBaseRotation() const {
    return getRep().defaultRB;
}
const Rotation& Constraint::ConstantOrientation::getDefaultFollowerRotation() const {
    return getRep().defaultRF;
}

    // ConstantOrientation bookkeeping //

bool Constraint::ConstantOrientation::isInstanceOf(const Constraint& s) {
    return ConstantOrientationRep::isA(s.getRep());
}
const Constraint::ConstantOrientation& Constraint::ConstantOrientation::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const ConstantOrientation&>(s);
}
Constraint::ConstantOrientation& Constraint::ConstantOrientation::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<ConstantOrientation&>(s);
}
const Constraint::ConstantOrientation::ConstantOrientationRep& Constraint::ConstantOrientation::getRep() const {
    return dynamic_cast<const ConstantOrientationRep&>(*rep);
}

Constraint::ConstantOrientation::ConstantOrientationRep& Constraint::ConstantOrientation::updRep() {
    return dynamic_cast<ConstantOrientationRep&>(*rep);
}

    // ConstantOrientationRep

    //TODO: no visualization yet



    //////////////////////
    // CONSTRAINT::WELD //
    //////////////////////

Constraint::Weld::Weld(MobilizedBody& body1, MobilizedBody& body2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Weld(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Weld(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new WeldRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B = updRep().addConstrainedBody(body1);
    updRep().F = updRep().addConstrainedBody(body2);
}

Constraint::Weld::Weld(MobilizedBody& body1, const Transform& frame1,
                       MobilizedBody& body2, const Transform& frame2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Weld(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Weld(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new WeldRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B = updRep().addConstrainedBody(body1);
    updRep().F = updRep().addConstrainedBody(body2);

    updRep().defaultFrameB = frame1;
    updRep().defaultFrameF = frame2;
}

Constraint::Weld& Constraint::Weld::setDefaultFrameOnBody1(const Transform& f1) {
    updRep().defaultFrameB = f1;
    return *this;
}

Constraint::Weld& Constraint::Weld::setDefaultFrameOnBody2(const Transform& f2) {
    updRep().defaultFrameF = f2;
    return *this;
}

MobilizedBodyIndex Constraint::Weld::getBody1MobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().B);
}
MobilizedBodyIndex Constraint::Weld::getBody2MobilizedBodyIndex() const {
    return getRep().getMobilizedBodyIndexOfConstrainedBody(getRep().F);
}
const Transform& Constraint::Weld::getDefaultFrameOnBody1() const {
    return getRep().defaultFrameB;
}
const Transform& Constraint::Weld::getDefaultFrameOnBody2() const {
    return getRep().defaultFrameF;
}


    // Weld bookkeeping //

bool Constraint::Weld::isInstanceOf(const Constraint& s) {
    return WeldRep::isA(s.getRep());
}
const Constraint::Weld& Constraint::Weld::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Weld&>(s);
}
Constraint::Weld& Constraint::Weld::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Weld&>(s);
}
const Constraint::Weld::WeldRep& Constraint::Weld::getRep() const {
    return dynamic_cast<const WeldRep&>(*rep);
}
Constraint::Weld::WeldRep& Constraint::Weld::updRep() {
    return dynamic_cast<WeldRep&>(*rep);
}

    // WeldRep

void Constraint::Weld::WeldRep::realizeTopologyVirtual(State& s) const { 
#ifdef USE_OLD_CONSTRAINTS
    SimbodyMatterSubsystemRep& matter = 
        const_cast<WeldRep*>(this)->updMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(F);
    const Vec3& station1 = defaultFrameB.T();
    const Vec3& station2 = defaultFrameF.T();
    const RigidBodyNode& rbn1 = matter.getRigidBodyNode(mobilizedBody1);
    const RigidBodyNode& rbn2 = matter.getRigidBodyNode(mobilizedBody2);
    const RBStation s1 (rbn1, station1);
    const RBStation s1x(rbn1, station1+defaultFrameB.x());
    const RBStation s1y(rbn1, station1+defaultFrameB.y());       
    const RBStation s1z(rbn1, station1+defaultFrameB.z());

    const RBStation s2 (rbn2, station2);
    const RBStation s2x(rbn2, station2+defaultFrameF.x());
    const RBStation s2y(rbn2, station2+defaultFrameF.y());       
    const RBStation s2z(rbn2, station2+defaultFrameF.z());

    // This is a "coincident station" constraint holding the frame
    // origins together (see above).
    matter.addOneDistanceConstraintEquation(s1x,s2,1.);
    matter.addOneDistanceConstraintEquation(s1y,s2,1.);
    matter.addOneDistanceConstraintEquation(s1z,s2,1.);

    // This is an "align axes" constraint. This has an unfortunate
    // symmetry when rotating 180 degrees about any axis.
    // This set of constraint equations is fine for *projection* but
    // not enough for *assembly*. You need to add another one to
    // eliminate the rotational symmetries when assembling from
    // far away.
    const Real d = std::sqrt(2.);
    matter.addOneDistanceConstraintEquation(s1y,s2z,d); // restrain x rot
    matter.addOneDistanceConstraintEquation(s1z,s2x,d); // restrain y rot
    matter.addOneDistanceConstraintEquation(s1x,s2y,d); // restrain z rot  
#endif
}

void Constraint::Weld::WeldRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the frames until we know the axis lengths to use, and we can't place
    // the geometry on the bodies until we know the body1 and body2 frame
    // placements, which might not be until Instance stage.
    if (stage == Stage::Instance && getAxisDisplayLength() != 0 && getMyMatterSubsystemRep().getShowDefaultGeometry()) {

        geom.push_back(DecorativeFrame(getAxisDisplayLength())
                                            .setColor(getFrameColor(0))
                                            .setLineThickness(2)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B))
                                            .setTransform(defaultFrameB));

        // Draw connector line back to body origin.
        if (defaultFrameB.T().norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultFrameB.T())
                             .setColor(getFrameColor(0))
                             .setLineThickness(2)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B)));
           

        geom.push_back(DecorativeFrame(0.67*getAxisDisplayLength())
                                            .setColor(getFrameColor(1))
                                            .setLineThickness(4)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(F))
                                            .setTransform(defaultFrameF));

        if (defaultFrameF.T().norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultFrameF.T())
                             .setColor(getFrameColor(1))
                             .setLineThickness(4)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(F)));
    }
}

    ////////////////////
    // CONSTRAINT REP //
    ////////////////////

/*virtual*/ Constraint::ConstraintRep::~ConstraintRep() {
    // NOTHING
}

void Constraint::ConstraintRep::realizeTopology(State& s) const
{
    // Calculate the relevant Subtree.
    mySubtree.clear();
    mySubtree.setSimbodyMatterSubsystem(getMyMatterSubsystem());
    for (ConstrainedBodyIndex b(0); b < (int)myConstrainedBodies.size(); ++b)
        mySubtree.addTerminalBody(myConstrainedBodies[b]);
    mySubtree.realizeTopology();

    realizeTopologyVirtual(s); // delegate to concrete constraint
}

void Constraint::ConstraintRep::invalidateTopologyCache() const {
    if (myMatterSubsystemRep)
        myMatterSubsystemRep->invalidateSubsystemTopologyCache();
}

bool Constraint::ConstraintRep::subsystemTopologyHasBeenRealized() const {
    return myMatterSubsystemRep && myMatterSubsystemRep->subsystemTopologyHasBeenRealized();
}

void Constraint::ConstraintRep::setMyMatterSubsystem
   (SimbodyMatterSubsystem& matter, ConstraintIndex id)
{
    assert(!isInSubsystem());
    myMatterSubsystemRep = &matter.updRep();
    myConstraintIndex = id;
}

const SimbodyMatterSubsystem& 
Constraint::ConstraintRep::getMyMatterSubsystem() const {
    return getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}

const MobilizedBody& 
Constraint::ConstraintRep::getConstrainedMobilizedBody(ConstrainedBodyIndex B) const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "Constrained bodies are not available until Topology stage has been realized.");
    return getMyMatterSubsystemRep().getMobilizedBody(myConstrainedBodies[B]);
}

const MobilizedBody& 
Constraint::ConstraintRep::getAncestorMobilizedBody() const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "The ancestor body is not available until Topology stage has been realized.");
    return getMyMatterSubsystemRep().getMobilizedBody(mySubtree.getAncestorMobilizedBodyIndex()); ;
}

// These are measured from and expressed in the ancestor (A) frame.
//TODO: should precalculate in State, return reference
Transform Constraint::ConstraintRep::getBodyTransform(const State& s, const SBPositionCache& pc, ConstrainedBodyIndex B) const { // X_AB
    const Transform& X_GB = getMyMatterSubsystemRep().getBodyTransform(s, pc, myConstrainedBodies[B]);
    const Transform& X_GA = getMyMatterSubsystemRep().getBodyTransform(s, pc, mySubtree.getAncestorMobilizedBodyIndex());
    return ~X_GA*X_GB;
}

SpatialVec Constraint::ConstraintRep::getBodyVelocity(const State& s, const SBVelocityCache& vc, ConstrainedBodyIndex B) const { // V_AB
    const Transform&  X_GB = getMyMatterSubsystemRep().getBodyTransform(s, myConstrainedBodies[B]);
    const Transform&  X_GA = getMyMatterSubsystemRep().getBodyTransform(s, mySubtree.getAncestorMobilizedBodyIndex());
    const SpatialVec& V_GB = getMyMatterSubsystemRep().getBodyVelocity(s, vc, myConstrainedBodies[B]);
    const SpatialVec& V_GA = getMyMatterSubsystemRep().getBodyVelocity(s, vc, mySubtree.getAncestorMobilizedBodyIndex());
    const Vec3 p_AB_G     = X_GB.T() - X_GA.T();
    const Vec3 p_AB_G_dot = V_GB[1]  - V_GA[1];        // d/dt p taken in G

    const Vec3 w_AB_G = V_GB[0] - V_GA[0];             // relative angular velocity of B in A, exp. in G

    // To get d/dt p taken in A, get derivative in G and remove the contribution generated by
    // A's velocity in G.
    const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // time deriv of p in A, exp in G

    return ~X_GA.R() * SpatialVec(w_AB_G, v_AB_G);     // re-express in A
}

SpatialVec Constraint::ConstraintRep::getBodyAcceleration(const State& s, const SBAccelerationCache& ac, ConstrainedBodyIndex B) const { // A_AB
    const Vec3&       p_GB = getMyMatterSubsystemRep().getBodyTransform(s, myConstrainedBodies[B]).T();
    const Transform&  X_GA = getMyMatterSubsystemRep().getBodyTransform(s, mySubtree.getAncestorMobilizedBodyIndex());
    const Vec3&       p_GA = X_GA.T();
    const SpatialVec& V_GB = getMyMatterSubsystemRep().getBodyVelocity(s, myConstrainedBodies[B]);
    const SpatialVec& V_GA = getMyMatterSubsystemRep().getBodyVelocity(s, mySubtree.getAncestorMobilizedBodyIndex());
    const SpatialVec& A_GB = getMyMatterSubsystemRep().getBodyAcceleration(s, ac, myConstrainedBodies[B]);
    const SpatialVec& A_GA = getMyMatterSubsystemRep().getBodyAcceleration(s, ac, mySubtree.getAncestorMobilizedBodyIndex());
    const Vec3&       w_GA = V_GA[0];
    const Vec3&       w_GB = V_GB[0];
    const Vec3&       b_GA = A_GA[0];
    const Vec3&       b_GB = A_GB[0];

    const Vec3 p_AB_G        = p_GB     - p_GA;
    const Vec3 p_AB_G_dot    = V_GB[1]  - V_GA[1];      // d/dt p taken in G
    const Vec3 p_AB_G_dotdot = A_GB[1]  - A_GA[1];      // d^2/dt^2 taken in G

    const Vec3 w_AB_G     = w_GB - w_GA;                // relative angular velocity of B in A, exp. in G
    const Vec3 v_AB_G     = p_AB_G_dot - w_GA % p_AB_G; // d/dt p taken in A, exp in G

    const Vec3 w_AB_G_dot = b_GB - b_GA;         // d/dt of w_AB_G taken in G
    const Vec3 v_AB_G_dot = p_AB_G_dotdot - (b_GA % p_AB_G + w_GA % p_AB_G_dot); // d/dt v_AB_G taken in G

    // We have the derivative in G; change it to derivative in A by adding in contribution caused
    // by motion of G in A, that is w_AG X w_AB_G. (Note that w_AG=-w_GA.)
    const Vec3 b_AB_G = w_AB_G_dot - w_GA % w_AB_G;
    const Vec3 a_AB_G = v_AB_G_dot - w_GA % v_AB_G; // taken in A, exp. in G

    return ~X_GA.R() * SpatialVec(b_AB_G, a_AB_G); // taken in A, expressed in A
}

// Find out how many holonomic (position), nonholonomic (velocity),
// and acceleration-only constraint equations are generated by this Constraint.
void Constraint::ConstraintRep::getNumConstraintEquations
   (const State& s, int& mp, int& mv, int& ma) const 
{
	const SBModelCache&   mc = getModelCache(s);
	const ConstraintIndex ix = myConstraintIndex;

	mp = mc.mHolonomicEquationsInUse[ix];
	mv = mc.mNonholonomicEquationsInUse[ix];
	ma = mc.mAccelerationOnlyEquationsInUse[ix];
}

// Find the slots in the QErr, UErr and UDotErr/Multiplier arrays allocated for the
// equations of this Constraint.
void Constraint::ConstraintRep::getConstraintEquationSlots
   (const State& s, int& holo0, int& nonholo0, int& accOnly0) const
{
	const SBModelCache&   mc = getModelCache(s);
	const ConstraintIndex ix = myConstraintIndex;

    holo0    = mc.holoErrSegment[ix].offset;
    nonholo0 = mc.nHolonomicConstraintEquationsInUse 
               + mc.nonholoErrSegment[ix].offset;
    accOnly0 = mc.nHolonomicConstraintEquationsInUse + mc.nNonholonomicConstraintEquationsInUse 
               + mc.accOnlyErrSegment[ix].offset;
}


// Given a state realized to Position stage, extract the position constraint errors
// corresponding to this Constraint. The 'mp' argument is for sanity checking -- it
// is an error if that isn't an exact match for the current number of holonomic
// constraint equations generated by this Constraint. We expect that perr points
// to an array of at least mp elements that we can write on.
void Constraint::ConstraintRep::getPositionErrors(const State& s, int mp, Real* perr) const {
	const SBModelCache& mc = getModelCache(s);

	assert(mp == mc.mHolonomicEquationsInUse[myConstraintIndex]);
	assert(mp == mc.holoErrSegment[myConstraintIndex].length);

	// Find the offset to our first qerr in the ModelCache.
	const int firstQErr = mc.holoErrSegment[myConstraintIndex].offset;

	// Get all qerr's for the subsystem.
	const Vector& qerr = getMyMatterSubsystemRep().getQErr(s);

	// Copy out the qerr's belonging to this constraint.
	for (int i=0; i < mp; ++i)
		perr[i] = qerr[firstQErr + i];
}

// Given a State realized to Velocity stage, extract the velocity constraint errors
// corresponding to this Constraint. This includes velocity constraints which were
// produced by differentiation of holonomic (position) constraints, and nonholonomic
// constraints which are introduced at the velocity level. The 'mpv' argument is
// for sanity checking -- it is an error if that isn't an exact match for the
// current number of holonomic+nonholonomic (mp+mv) constraint equations generated
// by this Constraint. We expect that pverr points to an array of at least mp+mv
// elements that we can write on.
void Constraint::ConstraintRep::getVelocityErrors(const State& s, int mpv, Real* pverr) const {
	const SBModelCache& mc = getModelCache(s);

	assert(mpv ==   mc.mHolonomicEquationsInUse   [myConstraintIndex] 
				  + mc.mNonholonomicEquationsInUse[myConstraintIndex]);
	assert(mpv ==  mc.holoErrSegment[myConstraintIndex].length
                 + mc.nonholoErrSegment[myConstraintIndex].length);

	// Get referente to all uerr's for the subsystem.
	const Vector& uerr = getMyMatterSubsystemRep().getUErr(s);

	// Find the offset to our first uerr in the ModelCache.
	const int firstHoloErr = mc.holoErrSegment[myConstraintIndex].offset;
    const int mHolo        = mc.holoErrSegment[myConstraintIndex].length;

    for (int i=0; i < mHolo; ++i)
        pverr[i] = uerr[firstHoloErr+i];

    const int firstNonholoErr = mc.nHolonomicConstraintEquationsInUse // total for whole subsystem
                                + mc.nonholoErrSegment[myConstraintIndex].offset;
    const int mNonholo        = mc.nonholoErrSegment[myConstraintIndex].length;

    for (int i=0; i < mNonholo; ++i)
        pverr[mHolo+i] = uerr[firstNonholoErr+i];
}

// Given a State realized to Acceleration stage, extract the accleration constraint errors
// corresponding to this Constraint. This includes acceleration constraints which were
// produced by twice differentiation of holonomic (position) constraints, and differentiation
// of nonholonomic (velocity) constraints, and acceleration-only constraints which are
// first introduced at the acceleration level. The 'mpva' argument is
// for sanity checking -- it is an error if that isn't an exact match for the
// current number of holonomic+nonholonomic+accelerationOnly (mp+mv+ma) constraint
// equations generated by this Constraint. We expect that pvaerr points to an array
// of at least mp+mv+ma elements that we can write on.
void Constraint::ConstraintRep::getAccelerationErrors(const State& s, int mpva, Real* pvaerr) const {
	const SBModelCache& mc = getModelCache(s);

	assert(mpva ==   mc.mHolonomicEquationsInUse       [myConstraintIndex] 
				   + mc.mNonholonomicEquationsInUse    [myConstraintIndex]
				   + mc.mAccelerationOnlyEquationsInUse[myConstraintIndex]);
	assert(mpva ==   mc.holoErrSegment[myConstraintIndex].length
                   + mc.nonholoErrSegment[myConstraintIndex].length
                   + mc.accOnlyErrSegment[myConstraintIndex].length);

	// Get referente to all udoterr's for the subsystem.
	const Vector& udoterr = getMyMatterSubsystemRep().getUDotErr(s);

	// Find the offset to our first uerr in the ModelCache.
	const int firstHoloErr = mc.holoErrSegment[myConstraintIndex].offset;
    const int mHolo        = mc.holoErrSegment[myConstraintIndex].length;

    for (int i=0; i < mHolo; ++i)
        pvaerr[i] = udoterr[firstHoloErr+i];

    const int firstNonholoErr = mc.nHolonomicConstraintEquationsInUse // total for whole subsystem
                                + mc.nonholoErrSegment[myConstraintIndex].offset;
    const int mNonholo        = mc.nonholoErrSegment[myConstraintIndex].length;

    for (int i=0; i < mNonholo; ++i)
        pvaerr[mHolo+i] = udoterr[firstNonholoErr+i];

    const int firstAccOnlyErr = mc.nHolonomicConstraintEquationsInUse+mc.nNonholonomicConstraintEquationsInUse // total for whole subsystem
                                + mc.accOnlyErrSegment[myConstraintIndex].offset;
    const int mAccOnly        = mc.accOnlyErrSegment[myConstraintIndex].length;

    for (int i=0; i < mAccOnly; ++i)
        pvaerr[mHolo+mNonholo+i] = udoterr[firstAccOnlyErr+i];
}

const SBModelCache& Constraint::ConstraintRep::getModelCache(const State& s) const {
    return getMyMatterSubsystemRep().getModelCache(s);
}
const SBPositionCache& Constraint::ConstraintRep::getPositionCache(const State& s) const {
    return getMyMatterSubsystemRep().getPositionCache(s);
}
const SBVelocityCache& Constraint::ConstraintRep::getVelocityCache(const State& s) const {
    return getMyMatterSubsystemRep().getVelocityCache(s);
}
const SBAccelerationCache& Constraint::ConstraintRep::getAccelerationCache(const State& s) const {
    return getMyMatterSubsystemRep().getAccelerationCache(s);
}

// Default implementations for ConstraintRep virtuals throw "unimplemented"
// exceptions. These shouldn't be called unless the concrete constraint has
// given a non-zero value for mp, mv, and/or ma which is a promise to 
// implement the associated methods.

    // These must be defined if there are any positin (holonomic) constraints defined.

void Constraint::ConstraintRep::
realizePositionErrorsVirtual(const State&, const SBPositionCache&, int mp,  Real* perr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "realizePositionErrors");
}

void Constraint::ConstraintRep::
realizePositionDotErrorsVirtual(const State&, const SBVelocityCache&, int mp,  Real* pverr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "realizePositionDotErrors");
}

void Constraint::ConstraintRep::
realizePositionDotDotErrorsVirtual(const State&, const SBAccelerationCache&, int mp,  Real* paerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "realizePositionDotDotErrors");
}


void Constraint::ConstraintRep::
applyPositionConstraintForcesVirtual
   (const State&, int mp, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "applyPositionConstraintForces");
}

    // These must be defined if there are any velocity (nonholonomic) constraints defined.

void Constraint::ConstraintRep::
realizeVelocityErrorsVirtual(const State&, const SBVelocityCache&, int mv,  Real* verr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "realizeVelocityErrors");
}


void Constraint::ConstraintRep::
realizeVelocityDotErrorsVirtual(const State&, const SBAccelerationCache&, int mv,  Real* vaerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "realizeVelocityDotErrors");
}


void Constraint::ConstraintRep::
applyVelocityConstraintForcesVirtual
   (const State&, int mv, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "applyVelocityConstraintForces");
}



// These must be defined if there are any acceleration-only constraints defined.
void Constraint::ConstraintRep::
realizeAccelerationErrorsVirtual(const State&, const SBAccelerationCache&, int ma,  Real* aerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "realizeAccelerationErrors");
}

void Constraint::ConstraintRep::
applyAccelerationConstraintForcesVirtual
   (const State&, int ma, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "applyAccelerationConstraintForces");
}


} // namespace SimTK

