/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-8 Stanford University and the Authors.         *
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

#include "ConstraintImpl.h"
#include "SimbodyMatterSubsystemRep.h"
#include "MobilizedBodyImpl.h"

#include <vector>
#include <algorithm>

#ifdef USE_OLD_CONSTRAINTS
    #include "LengthConstraints.h"
#endif

namespace SimTK {


    ////////////////
    // CONSTRAINT //
    ////////////////

void Constraint::disable(State& s) const {
    getImpl().setDisabled(s, true);
}
void Constraint::enable(State& s) const {
    getImpl().setDisabled(s, false);
}
bool Constraint::isDisabled(const State& s) const {
    return getImpl().isDisabled(s);
}
bool Constraint::isDisabledByDefault() const {
    return getImpl().isDisabledByDefault();
}
void Constraint::setDisabledByDefault(bool shouldBeDisabled) {
    updImpl().setDisabledByDefault(shouldBeDisabled);
}

const SimbodyMatterSubsystem& Constraint::getMatterSubsystem() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getMatterSubsystem() called on a Constraint that is not part of a subsystem.");
    return getImpl().getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}

ConstraintIndex Constraint::getConstraintIndex() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getConstraintIndex() called on a Constraint that is not part of a subsystem.");
    return getImpl().myConstraintIndex;
}

SimbodyMatterSubsystem& Constraint::updMatterSubsystem() {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "updMatterSubsystem() called on a Constraint that is not part of a subsystem.");
    return updImpl().updMyMatterSubsystemRep().updMySimbodyMatterSubsystemHandle();
}

bool Constraint::isInSubsystem() const {
    return getImpl().isInSubsystem();
}

bool Constraint::isInSameSubsystem(const MobilizedBody& body) const {
    return getImpl().isInSameSubsystem(body);
}

int Constraint::getNumConstrainedBodies() const {
    return getImpl().getNumConstrainedBodies();
}
int Constraint::getNumConstrainedMobilizers() const {
    return getImpl().getNumConstrainedMobilizers();
}

int Constraint::getNumConstrainedQ(const State& s) const {
    return getImpl().getNumConstrainedQ(s);
}
int Constraint::getNumConstrainedU(const State& s) const {
    return getImpl().getNumConstrainedU(s);
}

int Constraint::getNumConstrainedQ(const State& s, ConstrainedMobilizerIndex M) const {
    return getImpl().getNumConstrainedQ(s,M);
}
int Constraint::getNumConstrainedU(const State& s, ConstrainedMobilizerIndex M) const {
    return getImpl().getNumConstrainedU(s,M);
}

ConstrainedQIndex Constraint::getConstrainedQIndex
   (const State& s, ConstrainedMobilizerIndex M, MobilizerQIndex which) const {
    return getImpl().getConstrainedQIndex(s,M,which);
}
ConstrainedUIndex Constraint::getConstrainedUIndex
   (const State& s, ConstrainedMobilizerIndex M, MobilizerUIndex which) const {
    return getImpl().getConstrainedUIndex(s,M,which);
}

const MobilizedBody& Constraint::getMobilizedBodyFromConstrainedMobilizer(ConstrainedMobilizerIndex M) const {
    return getImpl().getMobilizedBodyFromConstrainedMobilizer(M);
}

const MobilizedBody& Constraint::getMobilizedBodyFromConstrainedBody(ConstrainedBodyIndex B) const {
    return getImpl().getMobilizedBodyFromConstrainedBody(B);
}
const MobilizedBody& Constraint::getAncestorMobilizedBody() const {
    return getImpl().getAncestorMobilizedBody();
}

const SimbodyMatterSubtree& Constraint::getSubtree() const {
    assert(getImpl().subsystemTopologyHasBeenRealized());
    return getImpl().mySubtree;
}

// Find out how many holonomic (position), nonholonomic (velocity),
// and acceleration-only constraint equations are generated by this Constraint.
void Constraint::getNumConstraintEquationsInUse(const State& s, int& mp, int& mv, int& ma) const {
	getImpl().getNumConstraintEquationsInUse(s,mp,mv,ma);
}

Vector Constraint::getPositionErrorsAsVector(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	Vector perr(mp);
	if (mp) getImpl().getPositionErrors(s, mp, &perr[0]);
	return perr;
}

Vector Constraint::getVelocityErrorsAsVector(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	Vector pverr(mp+mv);
	if (mp+mv) getImpl().getVelocityErrors(s, mp+mv, &pverr[0]);
	return pverr;
}

Vector Constraint::getAccelerationErrorsAsVector(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	Vector pvaerr(mp+mv+ma);
	if (mp+mv+ma) getImpl().getAccelerationErrors(s, mp+mv+ma, &pvaerr[0]);
	return pvaerr;
}

Vector Constraint::getMultipliersAsVector(const State& s) const {
    int mp,mv,ma;
    getNumConstraintEquationsInUse(s, mp, mv, ma);

    Vector mult(mp+mv+ma);
    if (mp+mv+ma) getImpl().getMultipliers(s, mp+mv+ma, &mult[0]);
    return mult;
}

Vector Constraint::calcPositionErrorFromQ(const State&, const Vector& q) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "Constraint", "calcPositionErrorFromQ");
}

Vector Constraint::calcVelocityErrorFromU(const State&, const Vector& q) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "Constraint", "calcVelocityErrorFromU");
}

Vector Constraint::calcAccelerationErrorFromUDot(const State&, const Vector& q) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "Constraint", "calcAccelerationErrorFromUDot");
}

Matrix Constraint::calcPositionConstraintMatrixP(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);

	Matrix P(mp, nu);
	if (mp && nu) {
		Vector  pverr0(mp), pverr(mp); // we're interested in the first mp of these
		State   tmp = s;      // don't change s

		matter.updU(tmp) = 0;	// first calculate the bias term -c(t,q)
		system.realize(tmp, Stage::Velocity);
		pverr0 = getVelocityErrorsAsVector(tmp)(0,mp);

		// Now calculate sensitivity of d(perr)/dt=Pu-c(t,q) to each u in turn.
		for (int j=0; j<nu; ++j) {
			matter.updU(tmp)[j] = 1;
		    system.realize(tmp, Stage::Velocity);
			pverr = getVelocityErrorsAsVector(tmp)(0,mp);
			matter.updU(tmp)[j] = 0;
			P(j) = pverr - pverr0;
		}
	}
	return P;
}

Matrix Constraint::calcPositionConstraintMatrixPt(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nb = matter.getNBodies();


	Matrix Pt(nu, mp);
    if (mp==0 || nu==0)
        return Pt;

    const ConstraintImpl& rep = getImpl();
    const int ncb = rep.getNumConstrainedBodies();
    const int ncu = rep.getNumConstrainedU(s);

    Vector              mobilityForces(ncu);
    Vector_<SpatialVec> bodyForcesInA(ncb); // might be zero of these

    Vector lambda(mp);
    lambda = 0;
    if (ncb == 0) {
        // Mobility forces only
        for (int i=0; i<mp; ++i) {
            lambda[i] = 1;
            mobilityForces = 0;
            rep.applyPositionConstraintForces(s, mp, &lambda[0], bodyForcesInA, mobilityForces);
            lambda[i] = 0;
            Pt(i) = 0;
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Pt(rep.getUIndexOfConstrainedU(s, cux), i) = mobilityForces[cux]; // unpack
        }
    } else {
        // There are some body forces
		Vector_<SpatialVec> bodyForcesInG(nb);
		bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

		// Calculate Pt*lambda with each lambda set to 1 in turn.
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
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Pt(rep.getUIndexOfConstrainedU(s, cux), i) += mobilityForces[cux]; // unpack
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
	getNumConstraintEquationsInUse(s, mp, mv, ma);

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
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nb = matter.getNBodies();


	Matrix Vt(nu, mv);
    if (mv==0 || nu==0)
        return Vt;

	const ConstraintImpl& rep = getImpl();
	const int ncb = rep.getNumConstrainedBodies();
    const int ncu = rep.getNumConstrainedU(s);

	Vector              mobilityForces(ncu);
	Vector_<SpatialVec> bodyForcesInA(ncb); // might be zero of these

	Vector lambda(mv);
    lambda = 0;
    if (ncb == 0) {
        // Mobility forces only
        for (int i=0; i<mv; ++i) {
		    lambda[i] = 1;
		    mobilityForces = 0;
		    rep.applyVelocityConstraintForces(s, mv, &lambda[0], bodyForcesInA, mobilityForces);
		    lambda[i] = 0;
            Vt(i) = 0;
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Vt(rep.getUIndexOfConstrainedU(s, cux), i) = mobilityForces[cux]; // unpack
        }
    } else {
        // There are some body forces
	    Vector_<SpatialVec> bodyForcesInG(nb);
	    bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

	    // Calculate Vt*lambda with each lambda set to 1 in turn.
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
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Vt(rep.getUIndexOfConstrainedU(s, cux), i) += mobilityForces[cux]; // unpack
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
	getNumConstraintEquationsInUse(s, mp, mv, ma);

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
	getNumConstraintEquationsInUse(s, mp, mv, ma);

	const SimbodyMatterSubsystem& matter = getMatterSubsystem();
	const System&                 system = matter.getSystem();

	const int nu = matter.getNU(s);
	const int nb = matter.getNBodies();


	Matrix At(nu, ma);
    if (ma==0 || nu==0)
        return At;

    const ConstraintImpl& rep = getImpl();
    const int ncb = rep.getNumConstrainedBodies();
    const int ncu = rep.getNumConstrainedU(s);

    Vector              mobilityForces(ncu);
    Vector_<SpatialVec> bodyForcesInA(ncb); // might be zero of these

    Vector lambda(ma);
    lambda = 0;
    
    if (ncb == 0) {
        // Mobility forces only
        for (int i=0; i<ma; ++i) {
            lambda[i] = 1;
            mobilityForces = 0;
            rep.applyAccelerationConstraintForces(s, ma, &lambda[0], bodyForcesInA, mobilityForces);
            lambda[i] = 0;
            At(i) = 0;
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                At(rep.getUIndexOfConstrainedU(s, cux), i) = mobilityForces[cux]; // unpack
        }
    } else {
        // There are some body forces
        Vector_<SpatialVec> bodyForcesInG(nb);
        bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

        // Calculate At*lambda with each lambda set to 1 in turn.
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
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                At(rep.getUIndexOfConstrainedU(s, cux), i) += mobilityForces[cux]; // unpack
        }
    }
	return At;
}

Matrix Constraint::calcPositionConstraintMatrixPQInverse(const State& s) const {
	int mp,mv,ma;
	getNumConstraintEquationsInUse(s, mp, mv, ma);

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
    getNumConstraintEquationsInUse(s, mp, mv, ma);
    assert(lambda.size() == mp+mv+ma);
    assert(lambda.hasContiguousData());

    getImpl().calcConstraintForcesFromMultipliers(s,mp,mv,ma,&lambda[0],bodyForcesInA,mobilityForces);
}


    /////////////////////
    // CONSTRAINT::ROD //
    /////////////////////

Constraint::Rod::Rod(MobilizedBody& body1, MobilizedBody& body2, Real defaultRodLength)
  : PIMPLDerivedHandleBase(new RodImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    //rep = new RodImpl(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);

    updImpl().defaultRodLength = defaultRodLength;
}

Constraint::Rod::Rod(MobilizedBody& body1, const Vec3& point1,
                     MobilizedBody& body2, const Vec3& point2, Real defaultRodLength)
  : PIMPLDerivedHandleBase(new RodImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    //rep = new RodImpl(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);

    updImpl().defaultPoint1 = point1;
    updImpl().defaultPoint2 = point2;
    updImpl().defaultRodLength = defaultRodLength;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody1(const Vec3& p1) {
    updImpl().defaultPoint1 = p1;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody2(const Vec3& p2) {
    updImpl().defaultPoint2 = p2;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultRodLength(Real length) {
    updImpl().defaultRodLength = length;
    return *this;
}


MobilizedBodyIndex Constraint::Rod::getBody1MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B1);
}
MobilizedBodyIndex Constraint::Rod::getBody2MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B1);
}
const Vec3& Constraint::Rod::getDefaultPointOnBody1() const {
    return getImpl().defaultPoint1;
}
const Vec3& Constraint::Rod::getDefaultPointOnBody2() const {
    return getImpl().defaultPoint2;
}
Real Constraint::Rod::getDefaultRodLength() const {
    return getImpl().defaultRodLength;
}

Real Constraint::Rod::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Real Constraint::Rod::getVelocityError(const State& s) const {
    Real pverr;
    getImpl().getVelocityErrors(s, 1, &pverr);
    return pverr;
}

Real Constraint::Rod::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::Rod::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}



    // RodImpl

void Constraint::Rod::RodImpl::realizeTopologyVirtual(State& s) const { 
#ifdef USE_OLD_CONSTRAINTS
    SimbodyMatterSubsystemRep& matter = 
        const_cast<RodImpl*>(this)->updMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B1);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(B2);
    const RigidBodyNode& rbn1 = matter.getRigidBodyNode(mobilizedBody1);
    const RigidBodyNode& rbn2 = matter.getRigidBodyNode(mobilizedBody2);
    const RBStation s1(rbn1, defaultPoint1);
    const RBStation s2(rbn2, defaultPoint2);
    matter.addOneDistanceConstraintEquation(s1,s2,defaultRodLength);
#endif
}

void Constraint::Rod::RodImpl::calcDecorativeGeometryAndAppendVirtual
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

        const Vec3 p_GP1 = getMobilizedBodyFromConstrainedBody(B1)
                              .findStationLocationInGround(s, defaultPoint1);
        const Vec3 p_GP2 = getMobilizedBodyFromConstrainedBody(B2)
                              .findStationLocationInGround(s, defaultPoint2);

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
  : PIMPLDerivedHandleBase(new PointInPlaneImpl())
{
    SimTK_ASSERT_ALWAYS(planeBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::PointInPlane(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(planeBody.isInSameSubsystem(followerBody),
        "Constraint::PointInPlane(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    //rep = new PointInPlaneRep(); rep->setMyHandle(*this);
    planeBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().planeBody    = updImpl().addConstrainedBody(planeBody);
    updImpl().followerBody = updImpl().addConstrainedBody(followerBody);
    updImpl().defaultPlaneNormal   = defPlaneNormal;
    updImpl().defaultPlaneHeight   = defPlaneHeight;
    updImpl().defaultFollowerPoint = defFollowerPoint;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultPlaneNormal(const UnitVec3& n) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultPlaneNormal = n;
    return *this;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultPlaneHeight(Real h) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultPlaneHeight = h;
    return *this;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultFollowerPoint(const Vec3& p) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultFollowerPoint = p;
    return *this;
}

MobilizedBodyIndex Constraint::PointInPlane::getPlaneMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().planeBody);
}
MobilizedBodyIndex Constraint::PointInPlane::getFollowerMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().followerBody);
}
const UnitVec3& Constraint::PointInPlane::getDefaultPlaneNormal() const {
    return getImpl().defaultPlaneNormal;
}
Real Constraint::PointInPlane::getDefaultPlaneHeight() const {
    return getImpl().defaultPlaneHeight;
}
const Vec3& Constraint::PointInPlane::getDefaultFollowerPoint() const {
    return getImpl().defaultFollowerPoint;
}

Constraint::PointInPlane& Constraint::PointInPlane::setPlaneDisplayHalfWidth(Real h) {
    updImpl().setPlaneDisplayHalfWidth(h);
    return *this;
}
Constraint::PointInPlane& Constraint::PointInPlane::setPointDisplayRadius(Real r) {
    updImpl().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::PointInPlane::getPlaneDisplayHalfWidth() const {
    return getImpl().getPlaneDisplayHalfWidth();
}

Real Constraint::PointInPlane::getPointDisplayRadius() const {
    return getImpl().getPointDisplayRadius();
}

Real Constraint::PointInPlane::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Real Constraint::PointInPlane::getVelocityError(const State& s) const {
    Real pverr;
    getImpl().getVelocityErrors(s, 1, &pverr);
    return pverr;
}

Real Constraint::PointInPlane::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::PointInPlane::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}

    // PointInPlaneImpl

void Constraint::PointInPlane::PointInPlaneImpl::calcDecorativeGeometryAndAppendVirtual
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
  : PIMPLDerivedHandleBase(new PointOnLineImpl())
{
    SimTK_ASSERT_ALWAYS(lineBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::PointOnLine(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(lineBody.isInSameSubsystem(followerBody),
        "Constraint::PointOnLine(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    //rep = new PointOnLineRep(); rep->setMyHandle(*this);
    lineBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().lineBody     = updImpl().addConstrainedBody(lineBody);
    updImpl().followerBody = updImpl().addConstrainedBody(followerBody);
    updImpl().defaultLineDirection = defLineDirection;
    updImpl().defaultPointOnLine   = defPointOnLine;
    updImpl().defaultFollowerPoint = defFollowerPoint;
}

Constraint::PointOnLine& Constraint::PointOnLine::setDefaultLineDirection(const UnitVec3& z) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultLineDirection = z;
    return *this;
}

Constraint::PointOnLine& Constraint::PointOnLine::setDefaultPointOnLine(const Vec3& P) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultPointOnLine = P;
    return *this;
}

Constraint::PointOnLine& Constraint::PointOnLine::setDefaultFollowerPoint(const Vec3& S) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultFollowerPoint = S;
    return *this;
}

MobilizedBodyIndex Constraint::PointOnLine::getLineMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().lineBody);
}
MobilizedBodyIndex Constraint::PointOnLine::getFollowerMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().followerBody);
}
const UnitVec3& Constraint::PointOnLine::getDefaultLineDirection() const {
    return getImpl().defaultLineDirection;
}
const Vec3& Constraint::PointOnLine::getDefaultPointOnLine() const {
    return getImpl().defaultPointOnLine;
}
const Vec3& Constraint::PointOnLine::getDefaultFollowerPoint() const {
    return getImpl().defaultFollowerPoint;
}

Constraint::PointOnLine& Constraint::PointOnLine::setLineDisplayHalfLength(Real h) {
    updImpl().setLineDisplayHalfLength(h);
    return *this;
}
Constraint::PointOnLine& Constraint::PointOnLine::setPointDisplayRadius(Real r) {
    updImpl().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::PointOnLine::getLineDisplayHalfLength() const {
    return getImpl().getLineDisplayHalfLength();
}

Real Constraint::PointOnLine::getPointDisplayRadius() const {
    return getImpl().getPointDisplayRadius();
}

Vec2 Constraint::PointOnLine::getPositionErrors(const State& s) const {
    Vec2 perr;
    getImpl().getPositionErrors(s, 2, &perr[0]);
    return perr;
}

Vec2 Constraint::PointOnLine::getVelocityErrors(const State& s) const {
    Vec2 pverr;
    getImpl().getVelocityErrors(s, 2, &pverr[0]);
    return pverr;
}

Vec2 Constraint::PointOnLine::getAccelerationErrors(const State& s) const {
    Vec2 pvaerr;
    getImpl().getAccelerationErrors(s, 2, &pvaerr[0]);
    return pvaerr;
}

Vec2 Constraint::PointOnLine::getMultipliers(const State& s) const {
    Vec2 mult;
    getImpl().getMultipliers(s, 2, &mult[0]);
    return mult;
}


    // PointOnLineImpl

void Constraint::PointOnLine::PointOnLineImpl::calcDecorativeGeometryAndAppendVirtual
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
  : PIMPLDerivedHandleBase(new ConstantAngleImpl())
{
    SimTK_ASSERT_ALWAYS(baseBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::ConstantAngle(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(baseBody.isInSameSubsystem(followerBody),
        "Constraint::ConstantAngle(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    //rep = new ConstantAngleRep(); rep->setMyHandle(*this);
    baseBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B = updImpl().addConstrainedBody(baseBody);
    updImpl().F = updImpl().addConstrainedBody(followerBody);
    updImpl().defaultAxisB = defaultAxisOnB;
    updImpl().defaultAxisF = defaultAxisOnF;
    updImpl().defaultAngle = angle;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setDefaultBaseAxis(const UnitVec3& a) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultAxisB = a;
    return *this;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setDefaultFollowerAxis(const UnitVec3& a) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultAxisF = a;
    return *this;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setDefaultAngle(Real t) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultAngle = t;
    return *this;
}

MobilizedBodyIndex Constraint::ConstantAngle::getBaseMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B);
}
MobilizedBodyIndex Constraint::ConstantAngle::getFollowerMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().F);
}
const UnitVec3& Constraint::ConstantAngle::getDefaultBaseAxis() const {
    return getImpl().defaultAxisB;
}
const UnitVec3& Constraint::ConstantAngle::getDefaultFollowerAxis() const {
    return getImpl().defaultAxisF;
}
Real Constraint::ConstantAngle::getDefaultAngle() const {
    return getImpl().defaultAngle;
}

Constraint::ConstantAngle& Constraint::ConstantAngle::setAxisDisplayLength(Real l) {
    updImpl().axisLength = l;
    return *this;
}
Constraint::ConstantAngle& Constraint::ConstantAngle::setAxisDisplayWidth(Real w) {
    updImpl().axisThickness = w;
    return *this;
}

Real Constraint::ConstantAngle::getAxisDisplayLength() const {
    return getImpl().axisLength;
}

Real Constraint::ConstantAngle::getAxisDisplayWidth() const {
    return getImpl().axisThickness;
}

Real Constraint::ConstantAngle::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Real Constraint::ConstantAngle::getVelocityError(const State& s) const {
    Real pverr;
    getImpl().getVelocityErrors(s, 1, &pverr);
    return pverr;
}

Real Constraint::ConstantAngle::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::ConstantAngle::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}

    // ConstantAngleImpl

void Constraint::ConstantAngle::ConstantAngleImpl::calcDecorativeGeometryAndAppendVirtual
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
  : PIMPLDerivedHandleBase(new BallImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Ball(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Ball(): both bodies to be connected must be in the same MatterSubsystem.");

    //rep = new BallRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);
}

Constraint::Ball::Ball(MobilizedBody& body1, const Vec3& point1,
                       MobilizedBody& body2, const Vec3& point2)
  : PIMPLDerivedHandleBase(new BallImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Ball(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Ball(): both bodies to be connected must be in the same MatterSubsystem.");

    //rep = new BallRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);

    updImpl().defaultPoint1 = point1;
    updImpl().defaultPoint2 = point2;
}

Constraint::Ball& Constraint::Ball::setDefaultPointOnBody1(const Vec3& p1) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultPoint1 = p1;
    return *this;
}

Constraint::Ball& Constraint::Ball::setDefaultPointOnBody2(const Vec3& p2) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultPoint2 = p2;
    return *this;
}

MobilizedBodyIndex Constraint::Ball::getBody1MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B1);
}
MobilizedBodyIndex Constraint::Ball::getBody2MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B2);
}
const Vec3& Constraint::Ball::getDefaultPointOnBody1() const {
    return getImpl().defaultPoint1;
}
const Vec3& Constraint::Ball::getDefaultPointOnBody2() const {
    return getImpl().defaultPoint2;
}

Constraint::Ball& Constraint::Ball::setDefaultRadius(Real r) {
    getImpl().invalidateTopologyCache();
    updImpl().setDefaultRadius(r);
    return *this;
}

Real Constraint::Ball::getDefaultRadius() const {
    return getImpl().getDefaultRadius();
}

Vec3 Constraint::Ball::getPositionErrors(const State& s) const {
    Vec3 perr;
    getImpl().getPositionErrors(s, 3, &perr[0]);
    return perr;
}

Vec3 Constraint::Ball::getVelocityErrors(const State& s) const {
    Vec3 pverr;
    getImpl().getVelocityErrors(s, 3, &pverr[0]);
    return pverr;
}

Vec3 Constraint::Ball::getAccelerationErrors(const State& s) const {
    Vec3 pvaerr;
    getImpl().getAccelerationErrors(s, 3, &pvaerr[0]);
    return pvaerr;
}

Vec3 Constraint::Ball::getMultipliers(const State& s) const {
    Vec3 mult;
    getImpl().getMultipliers(s, 3, &mult[0]);
    return mult;
}

    // BallImpl

void Constraint::Ball::BallImpl::realizeTopologyVirtual(State& s) const { 
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

void Constraint::Ball::BallImpl::calcDecorativeGeometryAndAppendVirtual
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

        // Draw connector line back to body origin.
        if (defaultPoint1.norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultPoint1)
                             .setColor(Gray)
                             .setLineThickness(2)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B1)));

        // On the outboard body draw an orange mesh sphere at the ball radius.
        geom.push_back(DecorativeSphere(getDefaultRadius())
                                            .setColor(Orange)
                                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                                            .setOpacity(0.5)
                                            .setResolution(0.5)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B2))
                                            .setTransform(X_B2));

        // Draw connector line back to body origin.
        if (defaultPoint2.norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultPoint2)
                             .setColor(Gray)
                             .setLineThickness(2)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B2)));
    }
}

    //////////////////////////////////////
    // CONSTRAINT::CONSTANT ORIENTATION //
    //////////////////////////////////////

Constraint::ConstantOrientation::ConstantOrientation
   (MobilizedBody& baseBody,     const Rotation& defaultFrameOnB,
    MobilizedBody& followerBody, const Rotation& defaultFrameOnF)
  : PIMPLDerivedHandleBase(new ConstantOrientationImpl())
{
    SimTK_ASSERT_ALWAYS(baseBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::ConstantOrientation(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(baseBody.isInSameSubsystem(followerBody),
        "Constraint::ConstantOrientation(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    //rep = new ConstantOrientationRep(); rep->setMyHandle(*this);
    baseBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B = updImpl().addConstrainedBody(baseBody);
    updImpl().F = updImpl().addConstrainedBody(followerBody);
    updImpl().defaultRB = defaultFrameOnB;
    updImpl().defaultRF = defaultFrameOnF;
}

Constraint::ConstantOrientation& Constraint::ConstantOrientation::setDefaultBaseRotation(const Rotation& R) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultRB = R;
    return *this;
}

Constraint::ConstantOrientation& Constraint::ConstantOrientation::setDefaultFollowerRotation(const Rotation& R) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultRF = R;
    return *this;
}


MobilizedBodyIndex Constraint::ConstantOrientation::getBaseMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B);
}
MobilizedBodyIndex Constraint::ConstantOrientation::getFollowerMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().F);
}
const Rotation& Constraint::ConstantOrientation::getDefaultBaseRotation() const {
    return getImpl().defaultRB;
}
const Rotation& Constraint::ConstantOrientation::getDefaultFollowerRotation() const {
    return getImpl().defaultRF;
}

Vec3 Constraint::ConstantOrientation::getPositionErrors(const State& s) const {
    Vec3 perr;
    getImpl().getPositionErrors(s, 3, &perr[0]);
    return perr;
}

Vec3 Constraint::ConstantOrientation::getVelocityErrors(const State& s) const {
    Vec3 pverr;
    getImpl().getVelocityErrors(s, 3, &pverr[0]);
    return pverr;
}

Vec3 Constraint::ConstantOrientation::getAccelerationErrors(const State& s) const {
    Vec3 pvaerr;
    getImpl().getAccelerationErrors(s, 3, &pvaerr[0]);
    return pvaerr;
}

Vec3 Constraint::ConstantOrientation::getMultipliers(const State& s) const {
    Vec3 mult;
    getImpl().getMultipliers(s, 3, &mult[0]);
    return mult;
}


    // ConstantOrientationImpl

    //TODO: no visualization yet



    //////////////////////
    // CONSTRAINT::WELD //
    //////////////////////

Constraint::Weld::Weld(MobilizedBody& body1, MobilizedBody& body2)
  : PIMPLDerivedHandleBase(new WeldImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Weld(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Weld(): both bodies to be connected must be in the same MatterSubsystem.");

    //rep = new WeldRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B = updImpl().addConstrainedBody(body1);
    updImpl().F = updImpl().addConstrainedBody(body2);
}

Constraint::Weld::Weld(MobilizedBody& body1, const Transform& frame1,
                       MobilizedBody& body2, const Transform& frame2)
  : PIMPLDerivedHandleBase(new WeldImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Weld(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Weld(): both bodies to be connected must be in the same MatterSubsystem.");

    //rep = new WeldRep(); rep->setMyHandle(*this);
    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B = updImpl().addConstrainedBody(body1);
    updImpl().F = updImpl().addConstrainedBody(body2);

    updImpl().defaultFrameB = frame1;
    updImpl().defaultFrameF = frame2;
}

Constraint::Weld& Constraint::Weld::setDefaultFrameOnBody1(const Transform& f1) {
    updImpl().defaultFrameB = f1;
    return *this;
}

Constraint::Weld& Constraint::Weld::setDefaultFrameOnBody2(const Transform& f2) {
    updImpl().defaultFrameF = f2;
    return *this;
}

MobilizedBodyIndex Constraint::Weld::getBody1MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B);
}
MobilizedBodyIndex Constraint::Weld::getBody2MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().F);
}
const Transform& Constraint::Weld::getDefaultFrameOnBody1() const {
    return getImpl().defaultFrameB;
}
const Transform& Constraint::Weld::getDefaultFrameOnBody2() const {
    return getImpl().defaultFrameF;
}

Vec6 Constraint::Weld::getPositionErrors(const State& s) const {
    Vec6 perr;
    getImpl().getPositionErrors(s, 6, &perr[0]);
    return perr;
}

Vec6 Constraint::Weld::getVelocityErrors(const State& s) const {
    Vec6 pverr;
    getImpl().getVelocityErrors(s, 6, &pverr[0]);
    return pverr;
}

Vec6 Constraint::Weld::getAccelerationErrors(const State& s) const {
    Vec6 pvaerr;
    getImpl().getAccelerationErrors(s, 6, &pvaerr[0]);
    return pvaerr;
}

Vec6 Constraint::Weld::getMultipliers(const State& s) const {
    Vec6 mult;
    getImpl().getMultipliers(s, 6, &mult[0]);
    return mult;
}

    // WeldImpl

void Constraint::Weld::WeldImpl::realizeTopologyVirtual(State& s) const { 
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

void Constraint::Weld::WeldImpl::calcDecorativeGeometryAndAppendVirtual
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


    ////////////////////////////
    // CONSTRAINT::NO SLIP 1D //
    ////////////////////////////

Constraint::NoSlip1D::NoSlip1D
   (MobilizedBody& caseBody, const Vec3& P_C, const UnitVec3& n_C,
    MobilizedBody& movingBody0, MobilizedBody& movingBody1)
  : PIMPLDerivedHandleBase(new NoSlip1DImpl())
{
    SimTK_ASSERT_ALWAYS(caseBody.isInSubsystem() && movingBody0.isInSubsystem()&& movingBody1.isInSubsystem(),
        "Constraint::NoSlip1D(): all three bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(caseBody.isInSameSubsystem(movingBody0) && caseBody.isInSameSubsystem(movingBody1),
        "Constraint::NoSlip1D(): all three bodies must be in the same SimbodyMatterSubsystem.");

    //rep = new NoSlip1DRep(); rep->setMyHandle(*this);
    caseBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().caseBody    = updImpl().addConstrainedBody(caseBody);
    updImpl().movingBody0 = updImpl().addConstrainedBody(movingBody0);
    updImpl().movingBody1 = updImpl().addConstrainedBody(movingBody1);
    updImpl().defaultNoSlipDirection = n_C;
    updImpl().defaultContactPoint    = P_C;
}

Constraint::NoSlip1D& Constraint::NoSlip1D::setDefaultDirection(const UnitVec3& n) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultNoSlipDirection = n;
    return *this;
}

Constraint::NoSlip1D& Constraint::NoSlip1D::setDefaultContactPoint(const Vec3& p) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultContactPoint = p;
    return *this;
}

MobilizedBodyIndex Constraint::NoSlip1D::getCaseMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().caseBody);
}
MobilizedBodyIndex Constraint::NoSlip1D::getMovingBodyMobilizedBodyIndex(int which) const {
    assert(which==0 || which==1);
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(
        which==0 ? getImpl().movingBody0 : getImpl().movingBody1);
}
const UnitVec3& Constraint::NoSlip1D::getDefaultDirection() const {
    return getImpl().defaultNoSlipDirection;
}
const Vec3& Constraint::NoSlip1D::getDefaultContactPoint() const {
    return getImpl().defaultContactPoint;
}

Constraint::NoSlip1D& Constraint::NoSlip1D::setDirectionDisplayLength(Real l) {
    updImpl().setDirectionDisplayLength(l);
    return *this;
}
Constraint::NoSlip1D& Constraint::NoSlip1D::setPointDisplayRadius(Real r) {
    updImpl().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::NoSlip1D::getDirectionDisplayLength() const {
    return getImpl().getDirectionDisplayLength();
}

Real Constraint::NoSlip1D::getPointDisplayRadius() const {
    return getImpl().getPointDisplayRadius();
}

Real Constraint::NoSlip1D::getVelocityError(const State& s) const {
    Real pverr;
    getImpl().getVelocityErrors(s, 1, &pverr);
    return pverr;
}

Real Constraint::NoSlip1D::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::NoSlip1D::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}

    // NoSlip1DImpl

void Constraint::NoSlip1D::NoSlip1DImpl::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the direction and contact
    // point location, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        // TODO: should be instance-stage data from State rather than topological data
        // This makes x axis point along no-slip direction, origin at contact point
        const Transform X_CP(Rotation(defaultNoSlipDirection,XAxis), defaultContactPoint);

        const MobilizedBodyIndex caseMBId = getMobilizedBodyIndexOfConstrainedBody(caseBody);

        if (directionLength > 0) {
            // On the case body, draw a gray line in the no-slip direction, starting at contact point.
            geom.push_back(DecorativeLine(Vec3(0), Vec3(directionLength,0,0))
                                                .setColor(Gray)
                                                .setBodyId(caseMBId)
                                                .setTransform(X_CP));
        }
        if (pointRadius > 0) {
            // On the follower body draw an orange mesh sphere at the point radius.
            geom.push_back(DecorativeSphere(pointRadius)
                                                .setColor(Orange)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setResolution(0.5)
                                                .setBodyId(caseMBId)
                                                .setTransform(X_CP));
        }
    }
}


    ////////////////////////////////
    // CONSTRAINT::CONSTANT SPEED //
    ////////////////////////////////

// This picks one of the mobilities from a multiple-mobility mobilizer.
Constraint::ConstantSpeed::ConstantSpeed
   (MobilizedBody& mobilizer, MobilizerUIndex whichU, Real defaultSpeed)
  : PIMPLDerivedHandleBase(new ConstantSpeedImpl())
{
    SimTK_ASSERT_ALWAYS(mobilizer.isInSubsystem(),
        "Constraint::ConstantSpeed(): the mobilizer must already be in a SimbodyMatterSubsystem.");

    //rep = new ConstantSpeedRep(); rep->setMyHandle(*this);
    mobilizer.updMatterSubsystem().adoptConstraint(*this);

    updImpl().theMobilizer = updImpl().addConstrainedMobilizer(mobilizer);
    updImpl().whichMobility = whichU;
    updImpl().prescribedSpeed = defaultSpeed;
}

// This is for mobilizers with only 1 mobility.
Constraint::ConstantSpeed::ConstantSpeed(MobilizedBody& mobilizer, Real defaultSpeed)
  : PIMPLDerivedHandleBase(new ConstantSpeedImpl())
{
    SimTK_ASSERT_ALWAYS(mobilizer.isInSubsystem(),
        "Constraint::ConstantSpeed(): the mobilizer must already be in a SimbodyMatterSubsystem.");

    //rep = new ConstantSpeedRep(); rep->setMyHandle(*this);
    mobilizer.updMatterSubsystem().adoptConstraint(*this);

    updImpl().theMobilizer = updImpl().addConstrainedMobilizer(mobilizer);
    updImpl().whichMobility = MobilizerUIndex(0);
    updImpl().prescribedSpeed = defaultSpeed;
}

MobilizedBodyIndex Constraint::ConstantSpeed::getMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedMobilizer(getImpl().theMobilizer);
}
MobilizerUIndex Constraint::ConstantSpeed::getWhichU() const {
    return getImpl().whichMobility;
}
Real Constraint::ConstantSpeed::getDefaultSpeed() const {
    return getImpl().prescribedSpeed;
}

Real Constraint::ConstantSpeed::getVelocityError(const State& s) const {
    Real pverr;
    getImpl().getVelocityErrors(s, 1, &pverr);
    return pverr;
}

Real Constraint::ConstantSpeed::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::ConstantSpeed::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}


    // ConstantSpeedImpl
    // nothing yet

    ////////////////////////
    // CONSTRAINT::CUSTOM //
    ////////////////////////

// We are given an Implementation object which is already holding a CustomImpl
// object for us. We'll first take away ownership of the CustomImpl, then
// make the CustomImpl take over ownership of the Implementation object.
Constraint::Custom::Custom(Constraint::Custom::Implementation* implementation)
  : PIMPLDerivedHandleBase(implementation ? implementation->updImpl().removeOwnershipOfCustomImpl()
                                          : 0)
{
    SimTK_ASSERT_ALWAYS(implementation,
        "Constraint::Custom::Custom(): Implementation pointer was NULL.");

    // Now store the Implementation pointer in our CustomImpl. The Implementation
    // object retains its original pointer to the CustomImpl object so it can
    // operate as a proxy for the CustomImpl. However the Custom handle now owns the
    // CustomImpl and the CustomImpl owns the Implementation.
    updImpl().takeOwnershipOfImplementation(implementation);

    updImpl().updMyMatterSubsystemRep().adoptConstraint(*this);
}

const Constraint::Custom::Implementation&
Constraint::Custom::getImplementation() const {
    return getImpl().getImplementation();
}

Constraint::Custom::Implementation&
Constraint::Custom::updImplementation() {
    return updImpl().updImplementation();
}

    // Constraint::CustomImpl

// The Implementation object should already contain a pointer to this CustomImpl object.
void Constraint::CustomImpl::takeOwnershipOfImplementation(Custom::Implementation* userImpl) {
    assert(!implementation); // you can only do this once!
    assert(userImpl);
    const Custom::ImplementationImpl& impImpl = userImpl->getImpl();
    assert(&impImpl.getCustomImpl() == this && !impImpl.isOwnerOfCustomImpl());
    implementation = userImpl;
}  

    ////////////////////////////////////////
    // CONSTRAINT::CUSTOM::IMPLEMENTATION //
    ////////////////////////////////////////

// Default constructor allocates a CustomImpl object and saves it in the ImplementationImpl object.
// When this gets passed to a Custom handle we'll turn over ownership of the CustomImpl object
// to the Custom handle.
Constraint::Custom::Implementation::Implementation(SimbodyMatterSubsystem& matter) 
  : PIMPLHandle<Implementation,ImplementationImpl>(new ImplementationImpl(new CustomImpl())) 
{
    // We don't know the ConstraintIndex yet since this hasn't been adopted by the MatterSubsystem.
    updImpl().updCustomImpl().setMyMatterSubsystem(matter, ConstraintIndex());
}

Constraint::Custom::Implementation::Implementation(SimbodyMatterSubsystem& matter, int mp, int mv, int ma) 
  : PIMPLHandle<Implementation,ImplementationImpl>(new ImplementationImpl(new CustomImpl(mp,mv,ma))) 
{
     // We don't know the ConstraintIndex yet since this hasn't been adopted by the MatterSubsystem.
   updImpl().updCustomImpl().setMyMatterSubsystem(matter, ConstraintIndex());
}

const SimbodyMatterSubsystem& Constraint::Custom::Implementation::getMatterSubsystem() const {
    return getImpl().getCustomImpl().getMyMatterSubsystem();
}


Constraint::Custom::Implementation*
Constraint::Custom::Implementation::clone() const {
    // This will not retain its connection to a CustomImpl class if it had one.
    return cloneVirtual();
}

void Constraint::Custom::Implementation::invalidateTopologyCache() const {
    getImpl().getCustomImpl().invalidateTopologyCache();
}

Constraint::Custom::Implementation& Constraint::Custom::Implementation::
setDefaultNumConstraintEquations(int mp, int mv, int ma) {
    updImpl().updCustomImpl().setDefaultNumConstraintEquations(mp,mv,ma);
    return *this;
}

Constraint::Custom::Implementation& Constraint::Custom::Implementation::
setDisabledByDefault(bool shouldBeDisabled) {
    updImpl().updCustomImpl().setDisabledByDefault(shouldBeDisabled);
    return *this;
}

ConstrainedBodyIndex Constraint::Custom::Implementation::
addConstrainedBody(const MobilizedBody& mb) {
    return updImpl().updCustomImpl().addConstrainedBody(mb);
}
ConstrainedMobilizerIndex Constraint::Custom::Implementation::
addConstrainedMobilizer(const MobilizedBody& mb) {
    return updImpl().updCustomImpl().addConstrainedMobilizer(mb);
}

MobilizedBodyIndex Constraint::Custom::Implementation::
getMobilizedBodyIndexOfConstrainedBody(ConstrainedBodyIndex c) const {
    return getImpl().getCustomImpl().getMobilizedBodyIndexOfConstrainedBody(c);
}
MobilizedBodyIndex Constraint::Custom::Implementation::
getMobilizedBodyIndexOfConstrainedMobilizer(ConstrainedMobilizerIndex c) const {
    return getImpl().getCustomImpl().getMobilizedBodyIndexOfConstrainedMobilizer(c);
}

Real Constraint::Custom::Implementation::
getOneQ(const State& s, ConstrainedMobilizerIndex cmx, MobilizerQIndex mqx) const {
    return getImpl().getCustomImpl().getOneQ(s,cmx,mqx);
}

Real Constraint::Custom::Implementation::
getOneU(const State& s, ConstrainedMobilizerIndex cmx, MobilizerUIndex mux) const {
    return getImpl().getCustomImpl().getOneU(s,cmx,mux);
}


Real Constraint::Custom::Implementation::
getOneQDot(const State& s, ConstrainedMobilizerIndex cmx, MobilizerQIndex mqx, bool realizingVelocity) const {
    return getImpl().getCustomImpl().getOneQDot(s,cmx,mqx,realizingVelocity);
}

Real Constraint::Custom::Implementation::
getOneQDotDot(const State& s, ConstrainedMobilizerIndex cmx, MobilizerQIndex mqx, bool realizingAcceleration) const {
    return getImpl().getCustomImpl().getOneQDotDot(s,cmx,mqx,realizingAcceleration);
}

Real Constraint::Custom::Implementation::
getOneUDot(const State& s, ConstrainedMobilizerIndex cmx, MobilizerUIndex mux, bool realizingAcceleration) const {
    return getImpl().getCustomImpl().getOneUDot(s,cmx,mux,realizingAcceleration);
}

// Apply a generalized (mobility) force to a particular mobility of the given constrained body B,
// adding it in to the appropriate slot of the mobilityForces vector.
void Constraint::Custom::Implementation::
addInOneMobilityForce(const State& s, ConstrainedMobilizerIndex M, MobilizerUIndex which,
                      Real f, Vector& mobilityForces) const 
{
    getImpl().getCustomImpl().addInOneMobilityForce(s,M,which,f,mobilityForces);
}

const Transform&  Constraint::Custom::Implementation::
getBodyTransform(const State& s, ConstrainedBodyIndex B, bool realizingPosition) const
{
    return getImpl().getCustomImpl().getBodyTransform(s,B,realizingPosition);
}

const SpatialVec& Constraint::Custom::Implementation::
getBodyVelocity(const State& s, ConstrainedBodyIndex B, bool realizingVelocity) const
{
    return getImpl().getCustomImpl().getBodyVelocity(s,B,realizingVelocity);
}

const SpatialVec& Constraint::Custom::Implementation::
getBodyAcceleration(const State& s, ConstrainedBodyIndex B, bool realizingAcceleration) const
{
    return getImpl().getCustomImpl().getBodyAcceleration(s,B,realizingAcceleration);
}

void Constraint::Custom::Implementation::
addInStationForce(const State& s, ConstrainedBodyIndex B, const Vec3& p_B, 
                       const Vec3& forceInA, Vector_<SpatialVec>& bodyForcesInA) const
{
    getImpl().getCustomImpl().addInStationForce(s,B,p_B,forceInA,bodyForcesInA);
}

void Constraint::Custom::Implementation::
addInBodyTorque(const State& s, ConstrainedBodyIndex B,
                     const Vec3& torqueInA, Vector_<SpatialVec>& bodyForcesInA) const
{
    getImpl().getCustomImpl().addInBodyTorque(s,B,torqueInA,bodyForcesInA);
}


// Default implementations for ConstraintImpl virtuals throw "unimplemented"
// exceptions. These shouldn't be called unless the concrete constraint has
// given a non-zero value for mp, mv, and/or ma which is a promise to 
// implement the associated methods.

    // These must be defined if there are any positin (holonomic) constraints defined.

void Constraint::Custom::Implementation::
realizePositionErrorsVirtual(const State&, int mp,  Real* perr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizePositionErrorsVirtual");
}

void Constraint::Custom::Implementation::
realizePositionDotErrorsVirtual(const State&, int mp,  Real* pverr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizePositionDotErrorsVirtual");
}

void Constraint::Custom::Implementation::
realizePositionDotDotErrorsVirtual(const State&, int mp,  Real* paerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizePositionDotDotErrorsVirtual");
}


void Constraint::Custom::Implementation::
applyPositionConstraintForcesVirtual
   (const State&, int mp, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "applyPositionConstraintForcesVirtual");
}

    // These must be defined if there are any velocity (nonholonomic) constraints defined.

void Constraint::Custom::Implementation::
realizeVelocityErrorsVirtual(const State&, int mv,  Real* verr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizeVelocityErrorsVirtual");
}


void Constraint::Custom::Implementation::
realizeVelocityDotErrorsVirtual(const State&, int mv,  Real* vaerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizeVelocityDotErrorsVirtual");
}


void Constraint::Custom::Implementation::
applyVelocityConstraintForcesVirtual
   (const State&, int mv, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "applyVelocityConstraintForcesVirtual");
}



// These must be defined if there are any acceleration-only constraints defined.
void Constraint::Custom::Implementation::
realizeAccelerationErrorsVirtual(const State&, int ma,  Real* aerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizeAccelerationErrorsVirtual");
}

void Constraint::Custom::Implementation::
applyAccelerationConstraintForcesVirtual
   (const State&, int ma, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "applyAccelerationConstraintForcesVirtual");
}

    ////////////////////////////////////
    // CONSTRAINT::COORDINATE COUPLER //
    ////////////////////////////////////

Constraint::CoordinateCoupler::CoordinateCoupler(SimbodyMatterSubsystem& matter, Function<1>* function, const std::vector<MobilizedBodyIndex>& coordBody, const std::vector<MobilizerQIndex>& coordIndex)
        : Custom(new CoordinateCouplerImpl(matter, function, coordBody, coordIndex)) {
}

Constraint::CoordinateCouplerImpl::CoordinateCouplerImpl(SimbodyMatterSubsystem& matter, Function<1>* function, const std::vector<MobilizedBodyIndex>& coordBody, const std::vector<MobilizerQIndex>& coordIndex)
        : Implementation(matter, 1, 0, 0), function(function), coordBodies(coordBody.size()), coordIndices(coordIndex), temp(coordBodies.size()), referenceCount(new int[1]) {
    assert(coordBodies.size() == coordIndices.size());
    assert(coordIndices.size() == function->getArgumentSize());
    assert(function->getMaxDerivativeOrder() >= 2);
    referenceCount[0] = 1;
    std::map<MobilizedBodyIndex,ConstrainedMobilizerIndex> bodyIndexMap;
    for (int i = 0; i < (int)coordBodies.size(); ++i) {
        if (bodyIndexMap.find(coordBody[i]) == bodyIndexMap.end())
            bodyIndexMap[coordBody[i]] = addConstrainedMobilizer(matter.getMobilizedBody(coordBody[i]));
        coordBodies[i] = bodyIndexMap[coordBody[i]];
    }
}

void Constraint::CoordinateCouplerImpl::realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQ(s, coordBodies[i], coordIndices[i]);
    perr[0] = function->calcValue(temp)[0];
}

void Constraint::CoordinateCouplerImpl::realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
    pverr[0] = 0.0;
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQ(s, coordBodies[i], coordIndices[i]);
    std::vector<int> components(1);
    for (int i = 0; i < temp.size(); ++i) {
        components[0] = i;
        pverr[0] += function->calcDerivative(components, temp)[0]*getOneQDot(s, coordBodies[i], coordIndices[i], true);
    }
}

void Constraint::CoordinateCouplerImpl::realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
    paerr[0] = 0.0;
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQ(s, coordBodies[i], coordIndices[i]);
    std::vector<int> components(2);
    // TODO this could be made faster by using symmetry if necessary
    for (int i = 0; i < temp.size(); ++i) {
        components[0] = i;
        Real qdoti = getOneQDot(s, coordBodies[i], coordIndices[i], true);
        for (int j = 0; j < temp.size(); ++j) {
            components[1] = j;
            paerr[0] += function->calcDerivative(components, temp)[0]*qdoti*getOneQDot(s, coordBodies[j], coordIndices[j], true);
        }
    }
    std::vector<int> component(1);
    const Vector& udot = s.updUDot();
    Vector qdotdot(s.getNQ());
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    SBStateDigest digest(s, matter.getRep(), Stage::Velocity);
    for (int i = 0; i < temp.size(); ++i) {
        component[0] = i;
        const MobilizedBody& body = matter.getMobilizedBody(getMobilizedBodyIndexOfConstrainedMobilizer(coordBodies[i]));
        const RigidBodyNode& node = body.getImpl().getMyRigidBodyNode();
        node.calcQDotDot(digest, udot, qdotdot);
        paerr[0] += function->calcDerivative(component, temp)[0]*body.getOneFromQPartition(s, coordIndices[i], qdotdot);
    }
}

void Constraint::CoordinateCouplerImpl::applyPositionConstraintForcesVirtual(const State& s, int mp, const Real* multipliers, Vector_<SpatialVec>& bodyForces, Vector& mobilityForces) const {
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQ(s, coordBodies[i], coordIndices[i]);
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    SBStateDigest digest(s, matter.getRep(), Stage::Velocity);
    bool useEuler = matter.getUseEulerAngles(s);
    std::vector<int> components(1);
    for (int i = 0; i < temp.size(); ++i) {
        components[0] = i;
        Real force = multipliers[0]*function->calcDerivative(components, temp)[0];
        const MobilizedBody& body = matter.getMobilizedBody(getMobilizedBodyIndexOfConstrainedMobilizer(coordBodies[i]));
        const RigidBodyNode& node = body.getImpl().getMyRigidBodyNode();
        Vector q = body.getQAsVector(s);
        Vector grad(body.getNumQ(s), 0.0);
        Vector forces(body.getNumU(s));
        grad[coordIndices[i]] = force;
        node.multiplyByQBlock(digest, useEuler, &q[0], true, &grad[0], &forces[0]);
        for (MobilizerUIndex index(0); index < forces.size(); index++)
            addInOneMobilityForce(s, coordBodies[i], index, forces[index], mobilityForces);
    }
}

    ///////////////////////////////
    // CONSTRAINT::SPEED COUPLER //
    ///////////////////////////////

Constraint::SpeedCoupler::SpeedCoupler(SimbodyMatterSubsystem& matter, Function<1>* function, const std::vector<MobilizedBodyIndex>& speedBody, const std::vector<MobilizerUIndex>& speedIndex)
        : Custom(new SpeedCouplerImpl(matter, function, speedBody, speedIndex, std::vector<MobilizedBodyIndex>(0), std::vector<MobilizerQIndex>(0))) {
}

Constraint::SpeedCoupler::SpeedCoupler(SimbodyMatterSubsystem& matter, Function<1>* function, const std::vector<MobilizedBodyIndex>& speedBody, const std::vector<MobilizerUIndex>& speedIndex,
        const std::vector<MobilizedBodyIndex>& coordBody, const std::vector<MobilizerQIndex>& coordIndex)
        : Custom(new SpeedCouplerImpl(matter, function, speedBody, speedIndex, coordBody, coordIndex)) {
}

Constraint::SpeedCouplerImpl::SpeedCouplerImpl(SimbodyMatterSubsystem& matter, Function<1>* function, const std::vector<MobilizedBodyIndex>& speedBody, const std::vector<MobilizerUIndex>& speedIndex,
        const std::vector<MobilizedBodyIndex>& coordBody, const std::vector<MobilizerQIndex>& coordIndex)
        : Implementation(matter, 0, 1, 0), function(function), speedBodies(speedBody.size()), speedIndices(speedIndex), coordBodies(coordBody), coordIndices(coordIndex),
        temp(speedBody.size()+coordBody.size()), referenceCount(new int[1]) {
    assert(speedBodies.size() == speedIndices.size());
    assert(coordBodies.size() == coordIndices.size());
    assert(temp.size() == function->getArgumentSize());
    assert(function->getMaxDerivativeOrder() >= 2);
    referenceCount[0] = 1;
    std::map<MobilizedBodyIndex,ConstrainedMobilizerIndex> bodyIndexMap;
    for (int i = 0; i < speedBodies.size(); ++i) {
        if (bodyIndexMap.find(speedBody[i]) == bodyIndexMap.end())
            bodyIndexMap[speedBody[i]] = addConstrainedMobilizer(matter.getMobilizedBody(speedBody[i]));
        speedBodies[i] = bodyIndexMap[speedBody[i]];
    }
}

void Constraint::SpeedCouplerImpl::realizeVelocityErrorsVirtual(const State& s, int mv,  Real* verr) const {
    findArguments(s);
    verr[0] = function->calcValue(temp)[0];
}

void Constraint::SpeedCouplerImpl::realizeVelocityDotErrorsVirtual(const State& s, int mv,  Real* vaerr) const {
    vaerr[0] = 0.0;
    findArguments(s);
    std::vector<int> components(1);
    for (int i = 0; i < (int) speedBodies.size(); ++i) {
        components[0] = i;
        vaerr[0] += function->calcDerivative(components, temp)[0]*getOneUDot(s, speedBodies[i], speedIndices[i], true);
    }
    for (int i = 0; i < (int) coordBodies.size(); ++i) {
        components[0] = speedBodies.size()+i;
        vaerr[0] += function->calcDerivative(components, temp)[0]*getMatterSubsystem().getMobilizedBody(coordBodies[i]).getOneQDot(s, coordIndices[i]);
    }
}

void Constraint::SpeedCouplerImpl::applyVelocityConstraintForcesVirtual(const State& s, int mv, const Real* multipliers, Vector_<SpatialVec>& bodyForces, Vector& mobilityForces) const {
    findArguments(s);
    std::vector<int> components(1);
    for (int i = 0; i < (int) speedBodies.size(); ++i) {
        components[0] = i;
        Real force = multipliers[0]*function->calcDerivative(components, temp)[0];
        addInOneMobilityForce(s, speedBodies[i], speedIndices[i], force, mobilityForces);
    }
}


    /////////////////////
    // CONSTRAINT IMPL //
    /////////////////////

void ConstraintImpl::realizeTopology(State& s) const
{
    // Calculate the relevant Subtree. There might not be any Constrained Bodies here
    // but we want to make sure we have a properly initialized empty Subtree in that case.
    mySubtree.clear();
    mySubtree.setSimbodyMatterSubsystem(getMyMatterSubsystem());
    for (ConstrainedBodyIndex b(0); b < myConstrainedBodies.size(); ++b)
        mySubtree.addTerminalBody(myConstrainedBodies[b]);
    mySubtree.realizeTopology();

    myAncestorBodyIsNotGround = myConstrainedBodies.size()
                                && mySubtree.getAncestorMobilizedBodyIndex() != GroundIndex;

    // If the Ancestor isn't Ground, reserve slots in the State cache
    // ancestor constrained body pools for each constrained body here
    // except the Ancestor (which may or may not be a constrained body).
    if (myAncestorBodyIsNotGround) {
        myPoolIndex.resize(myConstrainedBodies.size());
        for (ConstrainedBodyIndex b(0); b < myConstrainedBodies.size(); ++b) {
            if (myConstrainedBodies[b] == mySubtree.getAncestorMobilizedBodyIndex())
                myPoolIndex[b].invalidate();
            else myPoolIndex[b] = 
                getMyMatterSubsystemRep().allocateNextAncestorConstrainedBodyPoolSlot();
        }
    } else
        myPoolIndex.clear(); // ancestor is Ground; no need for pool entries

    realizeTopologyVirtual(s); // delegate to concrete constraint
}

// There are two main tasks here that can be performed now that we have values
// for the Model stage state variables:
// (1) Count up the number of holonomic, nonholonomic, and acceleration-only constraint
//     equations to be contributed by each Constraint, and assign corresponding slots
//     in constraint-equation ordered arrays, such as the State's constraint error arrays.
// (2) Above we assigned q's and u's to each mobilizer and stored the results in the
//     Model cache, now we can determine which of those q's and u's are involved in each
//     constraint. We need to collect up both the set of directly-constrained q's and u's
//     resulting from ConstrainedMobilizers, and indirectly-constrained ones arising 
//     from their effects on ConstrainedBodies. Together we call those "participating q's"
//     and "participating u's" (or "participating mobilities").
// The results of these computations goes in the Model cache.
void ConstraintImpl::realizeModel(State& s) const
{
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "ConstraintImpl::realizeModel() can't be called until after realizeToplogy().");

    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    const SBModelVars& modelVars  = matter.getModelVars(s);
    SBModelCache&      modelCache = matter.updModelCache(s);
    SBModelCache::PerConstraintModelInfo& cInfo =
        modelCache.updConstraintModelInfo(myConstraintIndex);

    cInfo.clear();
    cInfo.allocateConstrainedMobilizerModelInfo(getNumConstrainedMobilizers());

    if (isDisabled(s)) {
        cInfo.holoErrSegment    = Segment(0,modelCache.totalNHolonomicConstraintEquationsInUse);
        cInfo.nonholoErrSegment = Segment(0,modelCache.totalNNonholonomicConstraintEquationsInUse);
        cInfo.accOnlyErrSegment = Segment(0,modelCache.totalNAccelerationOnlyConstraintEquationsInUse);
        return;
    }

    // This constraint is not disabled

    // These are just the primary contraint equations, not their time derivatives.
    int mHolo, mNonholo, mAccOnly;
    calcNumConstraintEquationsInUse(s, mHolo, mNonholo, mAccOnly);

    // Must allocate space for the primary constraint equations and their time derivatives.
    //                                length         offset
    cInfo.holoErrSegment    = Segment(mHolo,    modelCache.totalNHolonomicConstraintEquationsInUse);
    cInfo.nonholoErrSegment = Segment(mNonholo, modelCache.totalNNonholonomicConstraintEquationsInUse);
    cInfo.accOnlyErrSegment = Segment(mAccOnly, modelCache.totalNAccelerationOnlyConstraintEquationsInUse);

    modelCache.totalNHolonomicConstraintEquationsInUse        += mHolo;
    modelCache.totalNNonholonomicConstraintEquationsInUse     += mNonholo;
    modelCache.totalNAccelerationOnlyConstraintEquationsInUse += mAccOnly;  

    // At this point we can find out how many q's and u's are associated with
    // each of the constrained mobilizers. We'll create packed arrays of q's and
    // u's ordered corresponding to the ConstrainedMobilizerIndices. We'll record
    // these in the ModelCache, by storing the ConstrainedQIndex and ConstrainedUIndex
    // of the lowest-numbered coordinate and mobility associated with each of
    // the ConstrainedMobilizers, along with the number of q's and u's.

    for (ConstrainedMobilizerIndex cmx(0); cmx < getNumConstrainedMobilizers(); ++cmx) {
        SBModelCache::PerConstrainedMobilizerModelInfo& mInfo = 
            cInfo.updConstrainedMobilizerModelInfo(cmx);

        const MobilizedBodyIndex mbx = getMobilizedBodyIndexOfConstrainedMobilizer(cmx);
        QIndex qix; int nq;
        UIndex uix; int nu;
        matter.findMobilizerQs(s,mbx,qix,nq);
        matter.findMobilizerUs(s,mbx,uix,nu);
        mInfo.nQInUse = nq;
        if (nq) {
            mInfo.firstConstrainedQIndex = cInfo.addConstrainedQ(qix);
            for (int i=1; i<nq; ++i) cInfo.addConstrainedQ(QIndex(qix+i));
        }
        mInfo.nUInUse = nu;
        if (nu) {
            mInfo.firstConstrainedUIndex = cInfo.addConstrainedU(uix);
            for (int i=1; i<nu; ++i) cInfo.addConstrainedU(UIndex(uix+i));
        }
    }

    // Now collect all the participating mobilities. This includes the constrained mobilities
    // as well as every q and u that can affect the constraint equations which involve 
    // constrained bodies. At the end we'll sort this list by subsystem QIndex/UIndex
    // and remove duplicates.
    cInfo.participatingQ = cInfo.constrainedQ;
    cInfo.participatingU = cInfo.constrainedU;

    const std::vector<MobilizedBodyIndex>& bodies = mySubtree.getAllBodies();
    for (int b=1; b<(int)bodies.size(); ++b) { // skip the Ancestor body 0
        QIndex qix; int nq;
        UIndex uix; int nu;
        matter.findMobilizerQs(s,bodies[b],qix,nq);
        matter.findMobilizerUs(s,bodies[b],uix,nu);
        for (int i=0; i<nq; ++i) cInfo.participatingQ.push_back(QIndex(qix+i));
        for (int i=0; i<nu; ++i) cInfo.participatingU.push_back(UIndex(uix+i));
    }

    std::sort(cInfo.participatingQ.begin(), cInfo.participatingQ.end());
    std::unique(cInfo.participatingQ.begin(), cInfo.participatingQ.end());

    realizeModelVirtual(s); // delegate to concrete constraint
}

void ConstraintImpl::realizeInstance(const State& s) const {
    if (isDisabled(s)) return;
    realizeInstanceVirtual(s); // nothing to do at the base class level
}
void ConstraintImpl::realizeTime(const State& s) const {
    if (isDisabled(s)) return;
    realizeTimeVirtual(s); // nothing to do in the base class
}
void ConstraintImpl::realizePosition(const State& s) const {
    if (isDisabled(s)) return;
    if (myAncestorBodyIsNotGround) {
        // Pre-calculate configuration information in the Ancestor frame
        SBPositionCache& pc = getMyMatterSubsystemRep().updPositionCache(s);
        for (ConstrainedBodyIndex b(0); b < getNumConstrainedBodies(); ++b) {
            const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[b];
            if (bx.isValid())
                pc.constrainedBodyConfigInAncestor[bx] = precalcConstrainedBodyTransformInAncestor(s,b);
        }
    }
    realizePositionVirtual(s); // delegate to concrete constraint
}
void ConstraintImpl::realizeVelocity(const State& s) const {
    if (isDisabled(s)) return;
    if (myAncestorBodyIsNotGround) {
        // Pre-calculate velocity information in the Ancestor frame
        SBVelocityCache& vc = getMyMatterSubsystemRep().updVelocityCache(s);
        for (ConstrainedBodyIndex b(0); b < getNumConstrainedBodies(); ++b) {
            const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[b];
            if (bx.isValid())
                vc.constrainedBodyVelocityInAncestor[bx] = precalcConstrainedBodyVelocityInAncestor(s,b);
        }
    }    
    realizeVelocityVirtual(s); // delegate to concrete constraint
}
void ConstraintImpl::realizeDynamics(const State& s) const {
    if (isDisabled(s)) return;
    realizeDynamicsVirtual(s); // nothing to do in the base class
}
void ConstraintImpl::realizeAcceleration(const State& s) const {
    if (isDisabled(s)) return;
    if (myAncestorBodyIsNotGround) {
        // Pre-calculate velocity information in the Ancestor frame
        SBAccelerationCache& ac = getMyMatterSubsystemRep().updAccelerationCache(s);
        for (ConstrainedBodyIndex b(0); b < getNumConstrainedBodies(); ++b) {
            const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[b];
            if (bx.isValid())
                ac.constrainedBodyAccelerationInAncestor[bx] = precalcConstrainedBodyAccelerationInAncestor(s,b);
        }
    }   
    realizeAccelerationVirtual(s); // delegate to concrete constraint
}
void ConstraintImpl::realizeReport(const State& s) const {
    if (isDisabled(s)) return;
    realizeReportVirtual(s); // nothing to do in the base class
}

void ConstraintImpl::invalidateTopologyCache() const {
    if (myMatterSubsystemRep)
        myMatterSubsystemRep->invalidateSubsystemTopologyCache();
}

bool ConstraintImpl::subsystemTopologyHasBeenRealized() const {
    return myMatterSubsystemRep && myMatterSubsystemRep->subsystemTopologyHasBeenRealized();
}

void ConstraintImpl::setMyMatterSubsystem
   (SimbodyMatterSubsystem& matter, ConstraintIndex id)
{
    // If this is already set it has to be to the same MatterSubsystem.
    assert(!myMatterSubsystemRep || myMatterSubsystemRep == &matter.getRep());
    myMatterSubsystemRep = &matter.updRep();
    myConstraintIndex = id;
}

const SimbodyMatterSubsystem& 
ConstraintImpl::getMyMatterSubsystem() const {
    return getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}

bool ConstraintImpl::isInSameSubsystem(const MobilizedBody& body) const {
    return isInSubsystem() && body.isInSubsystem() 
        && getMyMatterSubsystem().isSameSubsystem(body.getMatterSubsystem());
}

const MobilizedBody& 
ConstraintImpl::getMobilizedBodyFromConstrainedMobilizer(ConstrainedMobilizerIndex M) const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "Constrained mobilizers are not available until Topology stage has been realized.");
    return getMyMatterSubsystemRep().getMobilizedBody(myConstrainedMobilizers[M]);
}

const MobilizedBody& 
ConstraintImpl::getMobilizedBodyFromConstrainedBody(ConstrainedBodyIndex B) const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "Constrained bodies are not available until Topology stage has been realized.");
    return getMyMatterSubsystemRep().getMobilizedBody(myConstrainedBodies[B]);
}

const MobilizedBody& 
ConstraintImpl::getAncestorMobilizedBody() const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "The ancestor body is not available until Topology stage has been realized.");
    return getMyMatterSubsystemRep().getMobilizedBody(mySubtree.getAncestorMobilizedBodyIndex()); ;
}

Real ConstraintImpl::getOneQ
   (const State& s, ConstrainedMobilizerIndex cmx, MobilizerQIndex whichQ) const
{
    const QIndex qx = getQIndexOfConstrainedQ(s, getConstrainedQIndex(s, cmx, whichQ));
    return getMyMatterSubsystemRep().getQ(s)[qx];
}

Real ConstraintImpl::getOneU
   (const State& s, ConstrainedMobilizerIndex cmx, MobilizerUIndex whichU) const 
{
    const UIndex ux = getUIndexOfConstrainedU(s, getConstrainedUIndex(s, cmx, whichU));
    return getMyMatterSubsystemRep().getU(s)[ux];
}

Real ConstraintImpl::getOneQDot(const State& s, 
                ConstrainedMobilizerIndex cmx, MobilizerQIndex whichQ, bool realizing) const
{
    const QIndex qx = getQIndexOfConstrainedQ(s, getConstrainedQIndex(s, cmx, whichQ));
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return realizing ? matter.updQDot(s)[qx] : matter.getQDot(s)[qx];
}

Real ConstraintImpl::getOneUDot(const State& s,
                ConstrainedMobilizerIndex cmx, MobilizerUIndex whichU, bool realizing) const
{
    const UIndex ux = getUIndexOfConstrainedU(s, getConstrainedUIndex(s, cmx, whichU));
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return realizing ? matter.updUDot(s)[ux] : matter.getUDot(s)[ux];
}


Real ConstraintImpl::getOneQDotDot(const State& s, 
                ConstrainedMobilizerIndex cmx, MobilizerQIndex whichQ, bool realizing) const
{
    const QIndex qx = getQIndexOfConstrainedQ(s, getConstrainedQIndex(s, cmx, whichQ));
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return realizing ? matter.updQDotDot(s)[qx] : matter.getQDotDot(s)[qx];
}

// These are used to retrieve the indicated values from the State cache.
const Transform& ConstraintImpl::getBodyTransform
   (const State& s, ConstrainedBodyIndex B, bool realizingPosition) const 
{
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];

    if (!myAncestorBodyIsNotGround) 
        return matter.getBodyTransform(s,bodyB,realizingPosition); // X_GB

    static const Transform X_AA; // identity Transform
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();
    if (bodyB == ancestorA)
        return X_AA;

    const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[B];
    return matter.getPositionCache(s,realizingPosition)
                    .constrainedBodyConfigInAncestor[bx]; // X_AB
}
const SpatialVec& ConstraintImpl::getBodyVelocity
   (const State& s, ConstrainedBodyIndex B, bool realizingVelocity) const
{
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];

    if (!myAncestorBodyIsNotGround) 
        return matter.getBodyVelocity(s,bodyB,realizingVelocity); // V_GB

    static const SpatialVec V_AA(Vec3(0),Vec3(0)); // zero velocity
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();
    if (bodyB == ancestorA)
        return V_AA;

    const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[B];
    return matter.getVelocityCache(s,realizingVelocity)
                    .constrainedBodyVelocityInAncestor[bx]; // V_AB
}
const SpatialVec& ConstraintImpl::getBodyAcceleration
   (const State& s, ConstrainedBodyIndex B, bool realizingAcceleration) const
{
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];

    if (!myAncestorBodyIsNotGround) 
        return matter.getBodyAcceleration(s,bodyB,realizingAcceleration); // A_GB

    static const SpatialVec A_AA(Vec3(0),Vec3(0)); // zero velocity
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();
    if (bodyB == ancestorA)
        return A_AA;

    const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[B];
    return matter.getAccelerationCache(s,realizingAcceleration)
                    .constrainedBodyAccelerationInAncestor[bx]; // A_AB
}


// These are measured from and expressed in the ancestor (A) frame.
Transform ConstraintImpl::precalcConstrainedBodyTransformInAncestor(const State& s, ConstrainedBodyIndex B) const { // X_AB
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];
    const MobilizedBodyIndex         ancestorA = mySubtree.getAncestorMobilizedBodyIndex();

    const Transform& X_GB = matter.getBodyTransform(s, bodyB,     true);
    const Transform& X_GA = matter.getBodyTransform(s, ancestorA, true);
    return ~X_GA*X_GB;
}

SpatialVec ConstraintImpl::precalcConstrainedBodyVelocityInAncestor(const State& s, ConstrainedBodyIndex B) const { // V_AB
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];
    const MobilizedBodyIndex         ancestorA = mySubtree.getAncestorMobilizedBodyIndex();

    const Transform&  X_GB = matter.getBodyTransform(s, bodyB);
    const Transform&  X_GA = matter.getBodyTransform(s, ancestorA);
    const SpatialVec& V_GB = matter.getBodyVelocity(s, bodyB,     true);
    const SpatialVec& V_GA = matter.getBodyVelocity(s, ancestorA, true);
    const Vec3 p_AB_G     = X_GB.T() - X_GA.T();
    const Vec3 p_AB_G_dot = V_GB[1]  - V_GA[1];        // d/dt p taken in G

    const Vec3 w_AB_G = V_GB[0] - V_GA[0];             // relative angular velocity of B in A, exp. in G

    // To get d/dt p taken in A, get derivative in G and remove the contribution generated by
    // A's velocity in G.
    const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // time deriv of p in A, exp in G

    return ~X_GA.R() * SpatialVec(w_AB_G, v_AB_G);     // re-express in A
}

SpatialVec ConstraintImpl::precalcConstrainedBodyAccelerationInAncestor(const State& s, ConstrainedBodyIndex B) const { // A_AB
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];
    const MobilizedBodyIndex         ancestorA = mySubtree.getAncestorMobilizedBodyIndex();

    const Vec3&       p_GB = matter.getBodyTransform(s, bodyB).T();
    const Transform&  X_GA = matter.getBodyTransform(s, ancestorA);
    const Vec3&       p_GA = X_GA.T();
    const SpatialVec& V_GB = matter.getBodyVelocity(s, bodyB);
    const SpatialVec& V_GA = matter.getBodyVelocity(s, ancestorA);
    const SpatialVec& A_GB = matter.getBodyAcceleration(s, bodyB,     true);
    const SpatialVec& A_GA = matter.getBodyAcceleration(s, ancestorA, true);
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
// and acceleration-only constraint equations are generated by this Constraint
// as it is currently being modeled.
void ConstraintImpl::getNumConstraintEquationsInUse
   (const State& s, int& mHolo, int& mNonholo, int& mAccOnly) const 
{
    const SBModelCache::PerConstraintModelInfo& cInfo = 
        getModelCache(s).getConstraintModelInfo(myConstraintIndex);

	mHolo    = cInfo.holoErrSegment.length;
	mNonholo = cInfo.nonholoErrSegment.length;
	mAccOnly = cInfo.accOnlyErrSegment.length;
}

// Find the slots in the QErr, UErr and UDotErr/Multiplier arrays allocated for the
// equations of this Constraint.
void ConstraintImpl::getConstraintEquationSlots
   (const State& s, int& holo0, int& nonholo0, int& accOnly0) const
{
    const SBModelCache::PerConstraintModelInfo& cInfo = 
        getModelCache(s).getConstraintModelInfo(myConstraintIndex);

    const int mHolo    = cInfo.holoErrSegment.length;
	const int mNonholo = cInfo.nonholoErrSegment.length;

    holo0    =                    cInfo.holoErrSegment.offset;
    nonholo0 = mHolo            + cInfo.nonholoErrSegment.offset;
    accOnly0 = mHolo + mNonholo + cInfo.accOnlyErrSegment.offset;
}

void ConstraintImpl::setDisabled(State& s, bool shouldBeDisabled) const {
    getMyMatterSubsystemRep().setConstraintIsDisabled(s, myConstraintIndex, shouldBeDisabled);
}

bool ConstraintImpl::isDisabled(const State& s) const {
    return getMyMatterSubsystemRep().isConstraintDisabled(s, myConstraintIndex);
}

// Call this during construction phase to add a body to the topological structure of
// this Constraint. This body's mobilizer's mobilities are *not* part of the constraint; 
// mobilizers must be added separately.
ConstrainedBodyIndex ConstraintImpl::addConstrainedBody(const MobilizedBody& b) {
    assert(isInSameSubsystem(b));
    invalidateTopologyCache();

    const ConstrainedBodyIndex nextIx((int)myConstrainedBodies.size());

    // Add to the Mobilized->Constrained map and check for duplicates.
    std::pair<MobilizedBody2ConstrainedBodyMap::iterator, bool> result;
    result = myMobilizedBody2ConstrainedBodyMap.insert(
        MobilizedBody2ConstrainedBodyMap::value_type(b.getMobilizedBodyIndex(), nextIx));
    SimTK_ASSERT_ALWAYS(result.second,
        "addConstrainedBody(): a particular Constrained Body can be added only once per Constraint");

    // This is a new constrained body -- add it to the ConstrainedBody->MobilizedBody map too.
    myConstrainedBodies.push_back(b.getMobilizedBodyIndex());
    return nextIx;
}

// Call this during construction phase to add a mobilizer to the topological structure of
// this Constraint. All the coordinates q and mobilities u for this mobilizer are added also,
// but we don't know how many of those there will be until Stage::Model.
ConstrainedMobilizerIndex ConstraintImpl::addConstrainedMobilizer(const MobilizedBody& b) {
    assert(isInSameSubsystem(b));
    invalidateTopologyCache();

    const ConstrainedMobilizerIndex nextIx((int)myConstrainedMobilizers.size());

    // Add to the Mobilized->Constrained map and check for duplicates.
    std::pair<MobilizedBody2ConstrainedMobilizerMap::iterator, bool> result;
    result = myMobilizedBody2ConstrainedMobilizerMap.insert(
        MobilizedBody2ConstrainedMobilizerMap::value_type(b.getMobilizedBodyIndex(), nextIx));
    SimTK_ASSERT_ALWAYS(result.second,
        "addConstrainedMobilizer(): a particular Constrained Mobilizer can be added only once per Constraint");

    // This is a new constrained mobilizer -- add it to the ConstrainedMobilizer->MobilizedBody map too.
    myConstrainedMobilizers.push_back(b.getMobilizedBodyIndex());
    return nextIx;
}

QIndex ConstraintImpl::getQIndexOfConstrainedQ(const State& s, ConstrainedQIndex cqx) const {
    const SBModelCache&                         mc    = getModelCache(s);
    const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(myConstraintIndex);
    return cInfo.getQIndexFromConstrainedQ(cqx);
}

UIndex ConstraintImpl::getUIndexOfConstrainedU(const State& s, ConstrainedUIndex cqx) const {
    const SBModelCache&                         mc    = getModelCache(s);
    const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(myConstraintIndex);
    return cInfo.getUIndexFromConstrainedU(cqx);
}

int ConstraintImpl::getNumConstrainedQ(const State& s) const {
    return getModelCache(s).getConstraintModelInfo(myConstraintIndex).getNConstrainedQ();
}

int ConstraintImpl::getNumConstrainedQ
   (const State& s, ConstrainedMobilizerIndex M) const
{
    const SBModelCache::PerConstrainedMobilizerModelInfo& mInfo =
        getModelCache(s).getConstraintModelInfo(myConstraintIndex).getConstrainedMobilizerModelInfo(M);
    return mInfo.nQInUse; // same as corresponding MobilizedBody, or 0 if disabled
}

ConstrainedQIndex ConstraintImpl::getConstrainedQIndex
   (const State& s, ConstrainedMobilizerIndex M, MobilizerQIndex which) const 
{
    const int nq = getNumConstrainedQ(s,M);
    assert(0 <= which && which < nq);
    const SBModelCache::PerConstrainedMobilizerModelInfo& mInfo =
        getModelCache(s).getConstraintModelInfo(myConstraintIndex).getConstrainedMobilizerModelInfo(M);
    return ConstrainedQIndex(mInfo.firstConstrainedQIndex + which);
}       

int ConstraintImpl::getNumConstrainedU(const State& s) const {
    return getModelCache(s).getConstraintModelInfo(myConstraintIndex).getNConstrainedU();
}

int ConstraintImpl::getNumConstrainedU
   (const State& s, ConstrainedMobilizerIndex M) const
{
    const SBModelCache::PerConstrainedMobilizerModelInfo& mInfo =
        getModelCache(s).getConstraintModelInfo(myConstraintIndex).getConstrainedMobilizerModelInfo(M);
    return mInfo.nUInUse; // same as corresponding MobilizedBody, or 0 if disabled
}

ConstrainedUIndex ConstraintImpl::getConstrainedUIndex
   (const State& s, ConstrainedMobilizerIndex M, MobilizerUIndex which) const 
{
    const int nu = getNumConstrainedU(s,M);
    assert(0 <= which && which < nu);
    const SBModelCache::PerConstrainedMobilizerModelInfo& mInfo =
        getModelCache(s).getConstraintModelInfo(myConstraintIndex).getConstrainedMobilizerModelInfo(M);
    return ConstrainedUIndex(mInfo.firstConstrainedUIndex + which);
}   

// Given a state realized to Position stage, extract the position constraint errors
// corresponding to this Constraint. The 'mp' argument is for sanity checking -- it
// is an error if that isn't an exact match for the current number of holonomic
// constraint equations generated by this Constraint. We expect that perr points
// to an array of at least mp elements that we can write on.
void ConstraintImpl::getPositionErrors(const State& s, int mp, Real* perr) const {
    const SBModelCache::PerConstraintModelInfo& cInfo = 
        getModelCache(s).getConstraintModelInfo(myConstraintIndex);

	assert(mp == cInfo.holoErrSegment.length);

	// Find the offset to our first qerr in the ModelCache.
	const int firstQErr = cInfo.holoErrSegment.offset;

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
void ConstraintImpl::getVelocityErrors(const State& s, int mpv, Real* pverr) const {
    const SBModelCache& mc = getModelCache(s);
    const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(myConstraintIndex);

	assert(mpv ==  cInfo.holoErrSegment.length
                 + cInfo.nonholoErrSegment.length);

	// Get reference to all uerr's for the subsystem.
	const Vector& uerr = getMyMatterSubsystemRep().getUErr(s);

	// Find the offset to our first uerr in the ModelCache.
	const int firstHoloErr = cInfo.holoErrSegment.offset;
    const int mHolo        = cInfo.holoErrSegment.length;

    for (int i=0; i < mHolo; ++i)
        pverr[i] = uerr[firstHoloErr+i];

    const int firstNonholoErr = mc.totalNHolonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.nonholoErrSegment.offset;
    const int mNonholo        = cInfo.nonholoErrSegment.length;

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
void ConstraintImpl::getAccelerationErrors(const State& s, int mpva, Real* pvaerr) const {
    const SBModelCache& mc = getModelCache(s);
    const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(myConstraintIndex);

	assert(mpva ==   cInfo.holoErrSegment.length
                   + cInfo.nonholoErrSegment.length
                   + cInfo.accOnlyErrSegment.length);

	// Get reference to all udoterr's for the subsystem.
	const Vector& udoterr = getMyMatterSubsystemRep().getUDotErr(s);

	// Find the offset to our first uerr in the ModelCache.
	const int firstHoloErr = cInfo.holoErrSegment.offset;
    const int mHolo        = cInfo.holoErrSegment.length;

    for (int i=0; i < mHolo; ++i)
        pvaerr[i] = udoterr[firstHoloErr+i];

    const int firstNonholoErr = mc.totalNHolonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.nonholoErrSegment.offset;
    const int mNonholo        = cInfo.nonholoErrSegment.length;

    for (int i=0; i < mNonholo; ++i)
        pvaerr[mHolo+i] = udoterr[firstNonholoErr+i];

    const int firstAccOnlyErr = mc.totalNHolonomicConstraintEquationsInUse
                                + mc.totalNNonholonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.accOnlyErrSegment.offset;
    const int mAccOnly        = cInfo.accOnlyErrSegment.length;

    for (int i=0; i < mAccOnly; ++i)
        pvaerr[mHolo+mNonholo+i] = udoterr[firstAccOnlyErr+i];
}

// Given a State realized to Acceleration stage, extract the Lagrange multipliers
// corresponding to this Constraint. The 'mpva' argument is
// for sanity checking -- it is an error if that isn't an exact match for the
// current number of holonomic+nonholonomic+accelerationOnly (mp+mv+ma) constraint
// equations generated by this Constraint. We expect that lambda points to an array
// of at least mp+mv+ma elements that we can write on.
void ConstraintImpl::getMultipliers(const State& s, int mpva, Real* lambda) const {
    const SBModelCache& mc = getModelCache(s);
    const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(myConstraintIndex);

    assert(mpva ==   cInfo.holoErrSegment.length
                   + cInfo.nonholoErrSegment.length
                   + cInfo.accOnlyErrSegment.length);

    // Get reference to all multipliers for the subsystem.
    const Vector& multipliers = getMyMatterSubsystemRep().getMultipliers(s);

    // Find the offset to our first multiplier in the ModelCache.
    const int firstHoloErr = cInfo.holoErrSegment.offset;
    const int mHolo        = cInfo.holoErrSegment.length;

    for (int i=0; i < mHolo; ++i)
        lambda[i] = multipliers[firstHoloErr+i];

    const int firstNonholoErr = mc.totalNHolonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.nonholoErrSegment.offset;
    const int mNonholo        = cInfo.nonholoErrSegment.length;

    for (int i=0; i < mNonholo; ++i)
        lambda[mHolo+i] = multipliers[firstNonholoErr+i];

    const int firstAccOnlyErr = mc.totalNHolonomicConstraintEquationsInUse
                                + mc.totalNNonholonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.accOnlyErrSegment.offset;
    const int mAccOnly        = cInfo.accOnlyErrSegment.length;

    for (int i=0; i < mAccOnly; ++i)
        lambda[mHolo+mNonholo+i] = multipliers[firstAccOnlyErr+i];
}

const SBModelCache& ConstraintImpl::getModelCache(const State& s) const {
    return getMyMatterSubsystemRep().getModelCache(s);
}
const SBPositionCache& ConstraintImpl::getPositionCache(const State& s) const {
    return getMyMatterSubsystemRep().getPositionCache(s);
}
const SBVelocityCache& ConstraintImpl::getVelocityCache(const State& s) const {
    return getMyMatterSubsystemRep().getVelocityCache(s);
}
const SBAccelerationCache& ConstraintImpl::getAccelerationCache(const State& s) const {
    return getMyMatterSubsystemRep().getAccelerationCache(s);
}

// Default implementations for ConstraintImpl virtuals throw "unimplemented"
// exceptions. These shouldn't be called unless the concrete constraint has
// given a non-zero value for mp, mv, and/or ma which is a promise to 
// implement the associated methods.

    // These must be defined if there are any positin (holonomic) constraints defined.

void ConstraintImpl::
realizePositionErrorsVirtual(const State&, int mp,  Real* perr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "realizePositionErrors");
}

void ConstraintImpl::
realizePositionDotErrorsVirtual(const State&, int mp,  Real* pverr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "realizePositionDotErrors");
}

void ConstraintImpl::
realizePositionDotDotErrorsVirtual(const State&, int mp,  Real* paerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "realizePositionDotDotErrors");
}


void ConstraintImpl::
applyPositionConstraintForcesVirtual
   (const State&, int mp, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "applyPositionConstraintForces");
}

    // These must be defined if there are any velocity (nonholonomic) constraints defined.

void ConstraintImpl::
realizeVelocityErrorsVirtual(const State&, int mv,  Real* verr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "realizeVelocityErrors");
}


void ConstraintImpl::
realizeVelocityDotErrorsVirtual(const State&, int mv,  Real* vaerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "realizeVelocityDotErrors");
}


void ConstraintImpl::
applyVelocityConstraintForcesVirtual
   (const State&, int mv, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "applyVelocityConstraintForces");
}



// These must be defined if there are any acceleration-only constraints defined.
void ConstraintImpl::
realizeAccelerationErrorsVirtual(const State&, int ma,  Real* aerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "realizeAccelerationErrors");
}

void ConstraintImpl::
applyAccelerationConstraintForcesVirtual
   (const State&, int ma, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "applyAccelerationConstraintForces");
}


} // namespace SimTK

