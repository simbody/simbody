/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 *
 * Private implementation of Constraint, and its included subclasses which
 * represent the built-in constraint types.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintImpl.h"
#include "SimbodyMatterSubsystemRep.h"
#include "MobilizedBodyImpl.h"

#include <algorithm>

namespace SimTK {


//==============================================================================
//                               CONSTRAINT
//==============================================================================

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
    return getImpl().getMyConstraintIndex();
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

QIndex Constraint::getQIndexOfConstrainedQ(const State&      state,
                                           ConstrainedQIndex consQIndex) const {
    return getImpl().getQIndexOfConstrainedQ(state,consQIndex);
}

UIndex Constraint::getUIndexOfConstrainedU(const State&      state,
                                           ConstrainedUIndex consUIndex) const {
    return getImpl().getUIndexOfConstrainedU(state,consUIndex);
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
void Constraint::
getNumConstraintEquationsInUse(const State& s, int& mp, int& mv, int& ma) const 
{   getImpl().getNumConstraintEquationsInUse(s,mp,mv,ma); }

void Constraint::
getIndexOfMultipliersInUse(const State& s,
                           MultiplierIndex& px0, 
                           MultiplierIndex& vx0, 
                           MultiplierIndex& ax0) const
{   getImpl().getIndexOfMultipliersInUse(s,px0,vx0,ax0); }

void Constraint::setMyPartInConstraintSpaceVector(const State& state,
                                      const Vector& myPart,
                                      Vector& constraintSpace) const
{   getImpl().setMyPartInConstraintSpaceVector(state,myPart,constraintSpace); }

void Constraint::
getMyPartFromConstraintSpaceVector(const State& state,
                                   const Vector& constraintSpace,
                                   Vector& myPart) const
{   getImpl().getMyPartFromConstraintSpaceVector(state,constraintSpace,myPart);}

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

void Constraint::getConstraintForcesAsVectors
   (const State&         state,
    Vector_<SpatialVec>& bodyForcesInG,
    Vector&              mobilityForces) const 
{
    SimTK_ERRCHK1_ALWAYS(!isDisabled(state),
        "Constraint::getConstraintForcesAsVector()",
        "Constraint %d is currently disabled; you can't ask for its forces."
        " Use isDisabled() to check.", (int)getConstraintIndex());

    ArrayViewConst_<SpatialVec,ConstrainedBodyIndex>
        bodyF_G   = getImpl().getConstrainedBodyForcesInGFromState(state);
    ArrayViewConst_<Real,ConstrainedUIndex> 
        mobilityF = getImpl().getConstrainedMobilityForcesFromState(state);

    const int ncb = bodyF_G.size();
    bodyForcesInG.resize(ncb);
    for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx)
        bodyForcesInG[cbx] = bodyF_G[cbx];

    const int ncu = mobilityF.size();
    mobilityForces.resize(ncu);
    for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
        mobilityForces[cux] = mobilityF[cux];
}

// Multiply constraint forces (negated) by velocities to get power.
// (Remember that constraint forces are on the LHS so have the opposite
// sign from applied forces.)
Real Constraint::calcPower(const State& state) const {
    const ConstraintImpl& impl = getImpl();
    const SimbodyMatterSubsystemRep& matter = impl.getMyMatterSubsystemRep();

    ArrayViewConst_<SpatialVec,ConstrainedBodyIndex>
        bodyF_G   = impl.getConstrainedBodyForcesInGFromState(state);
    ArrayViewConst_<Real,ConstrainedUIndex> 
        mobilityF = impl.getConstrainedMobilityForcesFromState(state);
    
    Real power = 0;
    for (ConstrainedBodyIndex cbx(0); cbx < bodyF_G.size(); ++cbx) {
        const MobilizedBodyIndex mbx = 
            impl.getMobilizedBodyIndexOfConstrainedBody(cbx);
        const SpatialVec& V_GB = matter.getBodyVelocity(state, mbx);
        power -= ~bodyF_G[cbx] * V_GB;
    }
    const Vector& u = matter.getU(state);
    for (ConstrainedUIndex cux(0); cux < mobilityF.size(); ++cux) {
        const UIndex ux = impl.getUIndexOfConstrainedU(state, cux);
        power -= mobilityF[cux] * u[ux];
    }

    return power;
}


Vector Constraint::calcPositionErrorFromQ(const State&, const Vector& q) const {
    SimTK_THROW1(Exception::UnimplementedMethod, "Constraint::calcPositionErrorFromQ");
}

Vector Constraint::calcVelocityErrorFromU(const State&, const Vector& u) const {
    SimTK_THROW1(Exception::UnimplementedMethod, "Constraint::calcVelocityErrorFromU");
}

Vector Constraint::calcAccelerationErrorFromUDot(const State&, const Vector& udot) const {
    SimTK_THROW1(Exception::UnimplementedMethod, "Constraint::calcAccelerationErrorFromUDot");
}


// TODO: change to use operator form to avoid State copying kludge.
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

        matter.updU(tmp) = 0;   // first calculate the bias term -c(t,q)
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
    const int nb = matter.getNumBodies();


    Matrix Pt(nu, mp);
    if (mp==0 || nu==0)
        return Pt;

    const ConstraintImpl& rep = getImpl();
    const int ncb = rep.getNumConstrainedBodies();
    const int ncq = rep.getNumConstrainedQ(s);
    const int ncu = rep.getNumConstrainedU(s);

    // Any of these may be zero length.
    Array_<SpatialVec,ConstrainedBodyIndex> bodyForcesInA(ncb); 
    Array_<Real,ConstrainedQIndex>          qForces(ncq);
    Array_<Real,ConstrainedUIndex>          mobilityForces(ncu);

    Array_<Real> lambdap(mp, Real(0));

    if (ncb == 0) {
        // Mobility forces only
        for (int i=0; i<mp; ++i) {
            lambdap[i] = 1;
            qForces.fill(Real(0));
            rep.addInPositionConstraintForces(s, lambdap, 
                                              bodyForcesInA, qForces);           
            rep.convertQForcesToUForces(s, qForces, mobilityForces);  // fu = ~N fq
            lambdap[i] = 0;
            Pt(i) = 0;
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Pt(rep.getUIndexOfConstrainedU(s, cux), i) = mobilityForces[cux]; // unpack
        }
    } else {
        // There are some body forces
        Vector_<SpatialVec> bodyForcesInG(nb);
        bodyForcesInG.setToZero();

        // For converting those A-relative forces to G
        const Rotation& R_GA = rep.getAncestorMobilizedBody().getBodyRotation(s);

        // Calculate Pt*lambda with each lambda set to 1 in turn.
        for (int i=0; i<mp; ++i) {
            lambdap[i] = 1;
            bodyForcesInA.fill(SpatialVec(Vec3(0), Vec3(0)));
            qForces.fill(Real(0));
            rep.addInPositionConstraintForces(s, lambdap, 
                                              bodyForcesInA, qForces);
            rep.convertQForcesToUForces(s, qForces, mobilityForces);  // fu = ~N fq
            for (ConstrainedBodyIndex cb(0); cb < ncb; ++cb) {
                bodyForcesInG[rep.getMobilizedBodyIndexOfConstrainedBody(cb)] =
                    R_GA*bodyForcesInA[cb];
            }
            lambdap[i] = 0;

            rep.getMyMatterSubsystem().multiplyBySystemJacobianTranspose
                                                    (s,bodyForcesInG,Pt(i));
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
    const int nb = matter.getNumBodies();


    Matrix Vt(nu, mv);
    if (mv==0 || nu==0)
        return Vt;

    const ConstraintImpl& rep = getImpl();
    const int ncb = rep.getNumConstrainedBodies();
    const int ncu = rep.getNumConstrainedU(s);

        // Either of these may be zero length.
    Array_<SpatialVec,ConstrainedBodyIndex> bodyForcesInA(ncb); 
    Array_<Real,ConstrainedUIndex>          mobilityForces(ncu);
    Array_<Real> lambdav(mv, Real(0));

    if (ncb == 0) {
        // Mobility forces only
        for (int i=0; i<mv; ++i) {
            lambdav[i] = 1;
            mobilityForces.fill(Real(0));
            rep.addInVelocityConstraintForces(s, lambdav, 
                                              bodyForcesInA, mobilityForces);
            lambdav[i] = 0;
            Vt(i) = 0; // set column i to zero
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Vt(rep.getUIndexOfConstrainedU(s, cux), i) = 
                    mobilityForces[cux]; // unpack
        }
    } else {
        // There are some body forces
        Vector_<SpatialVec> bodyForcesInG(nb);
        bodyForcesInG = SpatialVec(Vec3(0), Vec3(0));

        // For converting those A-relative forces to G
        const Rotation& R_GA = 
            rep.getAncestorMobilizedBody().getBodyRotation(s);

        // Calculate Vt*lambda with each lambda set to 1 in turn.
        for (int i=0; i<mv; ++i) {
            lambdav[i] = 1;
            bodyForcesInA.fill(SpatialVec(Vec3(0), Vec3(0)));
            mobilityForces.fill(Real(0));
            rep.addInVelocityConstraintForces(s, lambdav, 
                                              bodyForcesInA, mobilityForces);
            for (ConstrainedBodyIndex cb(0); cb < ncb; ++cb) {
                bodyForcesInG[rep.getMobilizedBodyIndexOfConstrainedBody(cb)] =
                    R_GA*bodyForcesInA[cb];
            }
            lambdav[i] = 0;

            // Convert body forces F to generalized forces f=~J*F.
            // TODO: this should result in generalized forces only on the
            // participating mobilities - does it?
            rep.getMyMatterSubsystem().multiplyBySystemJacobianTranspose
                                                        (s,bodyForcesInG,Vt(i));

            // Now add in the generalized forces produced directly by the
            // Constraint, unpacking into the right global mobility slot.
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                Vt(rep.getUIndexOfConstrainedU(s,cux), i) 
                                                        += mobilityForces[cux];
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
    const int nb = matter.getNumBodies();


    Matrix At(nu, ma);
    if (ma==0 || nu==0)
        return At;

    const ConstraintImpl& rep = getImpl();
    const int ncb = rep.getNumConstrainedBodies();
    const int ncu = rep.getNumConstrainedU(s);

        // Either of these may be zero length.
    Array_<SpatialVec,ConstrainedBodyIndex> bodyForcesInA(ncb); 
    Array_<Real,ConstrainedUIndex>          mobilityForces(ncu);
    Array_<Real> lambdaa(ma, Real(0));
    
    if (ncb == 0) {
        // Mobility forces only
        for (int i=0; i<ma; ++i) {
            lambdaa[i] = 1;
            mobilityForces.fill(Real(0));
            rep.addInAccelerationConstraintForces(s, lambdaa, 
                                        bodyForcesInA, mobilityForces);
            lambdaa[i] = 0;
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
            lambdaa[i] = 1;
            bodyForcesInA.fill(SpatialVec(Vec3(0), Vec3(0)));
            mobilityForces.fill(Real(0));
            rep.addInAccelerationConstraintForces(s, lambdaa, 
                                        bodyForcesInA, mobilityForces);
            for (ConstrainedBodyIndex cb(0); cb < ncb; ++cb) {
                bodyForcesInG[rep.getMobilizedBodyIndexOfConstrainedBody(cb)] =
                    R_GA*bodyForcesInA[cb];
            }
            lambdaa[i] = 0;

            rep.getMyMatterSubsystem().multiplyBySystemJacobianTranspose
                                                        (s,bodyForcesInG,At(i));
            for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
                At(rep.getUIndexOfConstrainedU(s, cux), i) += mobilityForces[cux]; // unpack
        }
    }
    return At;
}

Matrix Constraint::calcPositionConstraintMatrixPNInv(const State& s) const {
    int mp,mv,ma;
    getNumConstraintEquationsInUse(s, mp, mv, ma);

    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    const System&                 system = matter.getSystem();

    const int nu = matter.getNU(s);
    const int nq = matter.getNQ(s);

    const Matrix P = calcPositionConstraintMatrixP(s);
    assert(P.nrow()==mp && P.ncol()==nu);

    Matrix PNInv(mp, nq); // = P*N^-1
    if (mp && nq) {
        // The routine below calculates qlikeRow = ulikeRow * N^-1 which is what
        // we need but it actually works with Vectors (transpose of RowVectors)
        // and at the moment they have to be contiguous so we'll have to copy.
        Vector uin(nu);
        Vector qout(nq);
        for (int i=0; i < mp; ++i) {
            uin = ~P[i];
            matter.multiplyByNInv(s, true, uin, qout);
            PNInv[i] = ~qout;
        }
    }
    return PNInv;
}

void Constraint::calcConstraintForcesFromMultipliers(
    const State&         s,
    const Vector&        lambda,
    Vector_<SpatialVec>& bodyForcesInA,
    Vector&              mobilityForces) const
{
    int mp, mv, ma;
    getNumConstraintEquationsInUse(s, mp, mv, ma);

    SimTK_APIARGCHECK2_ALWAYS(lambda.size() == mp+mv+ma,
        "Constraint", "calcConstraintForcesFromMultipliers()",
        "Expected %d multipliers but got %d.", mp+mv+ma, lambda.size());

    // We have to convert to and from Arrays since the underlying 
    // constraint methods use those for speed.

    Array_<Real> lambdap(mp), lambdav(mv), lambdaa(ma);
    for (int i=0; i<mp; ++i) lambdap[i] = lambda[i];
    for (int i=0; i<mv; ++i) lambdav[i] = lambda[mp+i];
    for (int i=0; i<ma; ++i) lambdaa[i] = lambda[mp+mv+i];

    Array_<SpatialVec,ConstrainedBodyIndex> tmpBodyForcesInA; 
    Array_<Real,      ConstrainedUIndex>    tmpMobilityForces;

    getImpl().calcConstraintForcesFromMultipliers
       (s,lambdap,lambdav,lambdaa,tmpBodyForcesInA,tmpMobilityForces);

    const int ncb = tmpBodyForcesInA.size();  // # constrained bodies
    const int ncu = tmpMobilityForces.size(); // # constrained u's

    bodyForcesInA.resize(ncb); mobilityForces.resize(ncu);

    for (ConstrainedBodyIndex i(0); i<ncb; ++i) 
        bodyForcesInA[i] = tmpBodyForcesInA[i];

    for (ConstrainedUIndex i(0); i<ncu; ++i)
        mobilityForces[i] = tmpMobilityForces[i];
}



//==============================================================================
//                             CONSTRAINT::ROD
//==============================================================================

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::Rod, Constraint::RodImpl, Constraint);

Constraint::Rod::Rod(MobilizedBody& body1, MobilizedBody& body2, Real defaultRodLength)
  : Constraint(new RodImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);

    updImpl().defaultRodLength = defaultRodLength;
}

Constraint::Rod::Rod(MobilizedBody& body1, const Vec3& point1,
                     MobilizedBody& body2, const Vec3& point2, Real defaultRodLength)
  : Constraint(new RodImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

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

void Constraint::Rod::RodImpl::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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
                                            .setColor(Gray)
                                            .setLineThickness(3)
                                            .setBodyId(GroundIndex));
        }
    }

}



//==============================================================================
//                         CONSTRAINT::POINT IN PLANE
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::PointInPlane, Constraint::PointInPlaneImpl, Constraint);

Constraint::PointInPlane::PointInPlane
   (MobilizedBody& planeBody,    const UnitVec3& defPlaneNormal, Real defPlaneHeight,
    MobilizedBody& followerBody, const Vec3&     defFollowerPoint)
  : Constraint(new PointInPlaneImpl())
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
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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



//==============================================================================
//                        CONSTRAINT::POINT ON LINE
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::PointOnLine, Constraint::PointOnLineImpl, Constraint);

Constraint::PointOnLine::PointOnLine
   (MobilizedBody& lineBody,     const UnitVec3& defLineDirection, const Vec3& defPointOnLine,
    MobilizedBody& followerBody, const Vec3&     defFollowerPoint)
  : Constraint(new PointOnLineImpl())
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
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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



//==============================================================================
//                         CONSTRAINT::CONSTANT ANGLE
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::ConstantAngle, Constraint::ConstantAngleImpl, Constraint);

Constraint::ConstantAngle::ConstantAngle
   (MobilizedBody& baseBody,     const UnitVec3& defaultAxisOnB,
    MobilizedBody& followerBody, const UnitVec3& defaultAxisOnF,
    Real angle)
  : Constraint(new ConstantAngleImpl())
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
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the normal, height, and follower
    // point location, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        //TODO
    }
}



//==============================================================================
//                              CONSTRAINT::BALL
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::Ball, Constraint::BallImpl, Constraint);

Constraint::Ball::Ball(MobilizedBody& body1, MobilizedBody& body2)
  : Constraint(new BallImpl())
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
  : Constraint(new BallImpl())
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

void Constraint::Ball::
setPointOnBody1(State& state, const Vec3& point_B1) const {
    getImpl().updBodyStations(state).first = point_B1;
}

void Constraint::Ball::
setPointOnBody2(State& state, const Vec3& point_B2) const {
    getImpl().updBodyStations(state).second = point_B2;
}

const Vec3& Constraint::Ball::
getPointOnBody1(const State& state) const {
    return getImpl().getBodyStations(state).first;
}
const Vec3& Constraint::Ball::
getPointOnBody2(const State& state) const {
    return getImpl().getBodyStations(state).second;
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

// The multipliers are the A-frame force vector as applied to body 2 at
// the user-defined station point on body 2. We want to return the
// force applied to body 1, expressed in body 1's frame. Note that this
// is the force applied to body 1 at the point coincident with body 2's
// station point -- if instead we returned the force at body 1's station
// point there would also be a small moment since that point in general
// differs from the contact point.
Vec3 Constraint::Ball::
getBallReactionForceOnBody1(const State& s) const {
    const BallImpl& impl = getImpl();
    Vec3 force2_A;
    impl.getMultipliers(s, 3, &force2_A[0]);
    const Rotation& R_AB1 = impl.getBodyRotationFromState(s, impl.B1);
    return -(~R_AB1*force2_A);
}

Vec3 Constraint::Ball::
getBallReactionForceOnBody2(const State& s) const{
    const BallImpl& impl = getImpl();
    Vec3 force2_A;
    impl.getMultipliers(s, 3, &force2_A[0]);
    const Rotation& R_AB2 = impl.getBodyRotationFromState(s, impl.B2);
    return ~R_AB2*force2_A;
}

    // BallImpl

// The default body stations may be overridden by setting instance variables
// in the state. We allocate the state resources here.
void Constraint::Ball::BallImpl::
realizeTopologyVirtual(State& state) const {
    stationsIx = getMyMatterSubsystemRep().
        allocateDiscreteVariable(state, Stage::Instance, 
            new Value< std::pair<Vec3,Vec3> >
               (std::make_pair(defaultPoint1,defaultPoint2)));
}

// Return the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Note that although
// these are used to define the position error, only the station on body 2
// is used to generate constraint forces; the point of body 1 that is 
// coincident with the body 2 point receives the equal and opposite force.
const std::pair<Vec3,Vec3>& Constraint::Ball::BallImpl::
getBodyStations(const State& state) const {
    return Value< std::pair<Vec3,Vec3> >::downcast
       (getMyMatterSubsystemRep().getDiscreteVariable(state,stationsIx));
}

// Return a writable reference into the Instance-stage state variable 
// containing the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Calling this
// method invalidates the Instance stage and above in the given state.
std::pair<Vec3,Vec3>& Constraint::Ball::BallImpl::
updBodyStations(State& state) const {
    return Value< std::pair<Vec3,Vec3> >::updDowncast
       (getMyMatterSubsystemRep().updDiscreteVariable(state,stationsIx));
}

void Constraint::Ball::BallImpl::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();

    // We can't generate the ball until we know the radius, and we can't place
    // the geometry on the body until we know the body1 and body2 point
    // placements on the bodies, which might not be until Instance stage.
    if (   stage == Stage::Instance 
        && getDefaultRadius() > 0 
        && matterRep.getShowDefaultGeometry()) 
    {
        const std::pair<Vec3,Vec3>& pts = getBodyStations(s);

        const Transform X_B1(Rotation(), pts.first); // should be point from State
        const Transform X_B2(Rotation(), pts.second);

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
        if (pts.first.norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), pts.first)
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
        if (pts.second.norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), pts.second)
                        .setColor(Gray)
                        .setLineThickness(2)
                        .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B2)));
    }
}



//==============================================================================
//                      CONSTRAINT::CONSTANT ORIENTATION
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::ConstantOrientation, Constraint::ConstantOrientationImpl, Constraint);

Constraint::ConstantOrientation::ConstantOrientation
   (MobilizedBody& baseBody,     const Rotation& defaultFrameOnB,
    MobilizedBody& followerBody, const Rotation& defaultFrameOnF)
  : Constraint(new ConstantOrientationImpl())
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



//==============================================================================
//                            CONSTRAINT::WELD
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::Weld, Constraint::WeldImpl, Constraint);

Constraint::Weld::Weld(MobilizedBody& body1, MobilizedBody& body2)
  : Constraint(new WeldImpl())
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
  : Constraint(new WeldImpl())
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

void Constraint::Weld::WeldImpl::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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
        if (defaultFrameB.p().norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultFrameB.p())
                             .setColor(getFrameColor(0))
                             .setLineThickness(2)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(B)));
           

        geom.push_back(DecorativeFrame(0.67*getAxisDisplayLength())
                                            .setColor(getFrameColor(1))
                                            .setLineThickness(4)
                                            .setBodyId(getMobilizedBodyIndexOfConstrainedBody(F))
                                            .setTransform(defaultFrameF));

        if (defaultFrameF.p().norm() >= SignificantReal)
            geom.push_back(DecorativeLine(Vec3(0), defaultFrameF.p())
                             .setColor(getFrameColor(1))
                             .setLineThickness(4)
                             .setBodyId(getMobilizedBodyIndexOfConstrainedBody(F)));
    }
}



//==============================================================================
//                            CONSTRAINT::NO SLIP 1D
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::NoSlip1D, Constraint::NoSlip1DImpl, Constraint);

Constraint::NoSlip1D::NoSlip1D
   (MobilizedBody& caseBody, const Vec3& P_C, const UnitVec3& n_C,
    MobilizedBody& movingBody0, MobilizedBody& movingBody1)
  : Constraint(new NoSlip1DImpl())
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
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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



//==============================================================================
//                         CONSTRAINT::CONSTANT SPEED
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::ConstantSpeed, Constraint::ConstantSpeedImpl, Constraint);

// This picks one of the mobilities from a multiple-mobility mobilizer.
Constraint::ConstantSpeed::ConstantSpeed
   (MobilizedBody& mobilizer, MobilizerUIndex whichU, Real defaultSpeed)
  : Constraint(new ConstantSpeedImpl())
{
    SimTK_ASSERT_ALWAYS(mobilizer.isInSubsystem(),
        "Constraint::ConstantSpeed(): the mobilizer must already be in a SimbodyMatterSubsystem.");

    mobilizer.updMatterSubsystem().adoptConstraint(*this);

    updImpl().theMobilizer = updImpl().addConstrainedMobilizer(mobilizer);
    updImpl().whichMobility = whichU;
    updImpl().prescribedSpeed = defaultSpeed;
}

// This is for mobilizers with only 1 mobility.
Constraint::ConstantSpeed::ConstantSpeed(MobilizedBody& mobilizer, Real defaultSpeed)
  : Constraint(new ConstantSpeedImpl())
{
    SimTK_ASSERT_ALWAYS(mobilizer.isInSubsystem(),
        "Constraint::ConstantSpeed(): the mobilizer must already be in a SimbodyMatterSubsystem.");

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



//==============================================================================
//                     CONSTRAINT::CONSTANT ACCELERATION
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS
   (Constraint::ConstantAcceleration, Constraint::ConstantAccelerationImpl, 
    Constraint);

// This picks one of the mobilities from a multiple-mobility mobilizer.
Constraint::ConstantAcceleration::ConstantAcceleration
   (MobilizedBody& mobilizer, MobilizerUIndex whichU, Real defaultAcceleration)
:   Constraint(new ConstantAccelerationImpl())
{
    SimTK_ASSERT_ALWAYS(mobilizer.isInSubsystem(),
    "Constraint::ConstantAcceleration(): the mobilizer must already be"
    " in a SimbodyMatterSubsystem.");

    mobilizer.updMatterSubsystem().adoptConstraint(*this);

    updImpl().theMobilizer = updImpl().addConstrainedMobilizer(mobilizer);
    updImpl().whichMobility = whichU;
    updImpl().defaultAcceleration = defaultAcceleration;
}

// This is for mobilizers with only 1 mobility.
Constraint::ConstantAcceleration::ConstantAcceleration
   (MobilizedBody& mobilizer, Real defaultAcceleration)
:   Constraint(new ConstantAccelerationImpl())
{
    SimTK_ASSERT_ALWAYS(mobilizer.isInSubsystem(),
    "Constraint::ConstantAcceleration(): the mobilizer must already be"
    " in a SimbodyMatterSubsystem.");

    mobilizer.updMatterSubsystem().adoptConstraint(*this);

    updImpl().theMobilizer = updImpl().addConstrainedMobilizer(mobilizer);
    updImpl().whichMobility = MobilizerUIndex(0);
    updImpl().defaultAcceleration = defaultAcceleration;
}

MobilizedBodyIndex Constraint::ConstantAcceleration::
getMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedMobilizer
                                            (getImpl().theMobilizer);
}
MobilizerUIndex Constraint::ConstantAcceleration::getWhichU() const {
    return getImpl().whichMobility;
}
Real Constraint::ConstantAcceleration::getDefaultAcceleration() const {
    return getImpl().defaultAcceleration;
}
Constraint::ConstantAcceleration& Constraint::ConstantAcceleration::
setDefaultAcceleration(Real udot) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultAcceleration = udot;
    return *this;
}

void Constraint::ConstantAcceleration::
setAcceleration(State& state, Real accel) const {
    getImpl().updAcceleration(state) = accel;
}

Real Constraint::ConstantAcceleration::
getAcceleration(const State& state) const {
    return getImpl().getAcceleration(state);
}

Real Constraint::ConstantAcceleration::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::ConstantAcceleration::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}

    // ConstantAccelerationImpl

// Allocate a state variable to hold the desired acceleration.
void Constraint::ConstantAccelerationImpl::
realizeTopologyVirtual(State& state) const {
    ConstantAccelerationImpl* mThis = // mutable momentarily
        const_cast<ConstantAccelerationImpl*>(this);
    mThis->accelIx = getMyMatterSubsystemRep().
        allocateDiscreteVariable(state, Stage::Acceleration, 
            new Value<Real>(defaultAcceleration));
}

Real Constraint::ConstantAccelerationImpl::
getAcceleration(const State& state) const {
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return Value<Real>::downcast(matter.getDiscreteVariable(state,accelIx));
}

Real& Constraint::ConstantAccelerationImpl::
updAcceleration(State& state) const {
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return Value<Real>::updDowncast(matter.updDiscreteVariable(state,accelIx));
}



//==============================================================================
//                            CONSTRAINT::CUSTOM
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS
   (Constraint::Custom, Constraint::CustomImpl, Constraint);

// We are given an Implementation object which is already holding a CustomImpl
// object for us. We'll first take away ownership of the CustomImpl, then
// make the CustomImpl take over ownership of the Implementation object.
Constraint::Custom::Custom(Constraint::Custom::Implementation* implementation)
:   Constraint(implementation 
                ? implementation->updImpl().removeOwnershipOfCustomImpl()
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
void Constraint::CustomImpl::
takeOwnershipOfImplementation(Custom::Implementation* userImpl) {
    assert(!implementation); // you can only do this once!
    assert(userImpl);
    const Custom::ImplementationImpl& impImpl = userImpl->getImpl();
    assert(&impImpl.getCustomImpl() == this && !impImpl.isOwnerOfCustomImpl());
    implementation = userImpl;
}  



//==============================================================================
//                     CONSTRAINT::CUSTOM::IMPLEMENTATION
//==============================================================================

// Default constructor allocates a CustomImpl object and saves it in the 
// ImplementationImpl object. When this gets passed to a Custom handle we'll 
// turn over ownership of the CustomImpl object to the Custom handle.
Constraint::Custom::Implementation::Implementation
   (SimbodyMatterSubsystem& matter) 
:   PIMPLHandle<Implementation,ImplementationImpl>
        (new ImplementationImpl(new CustomImpl())) 
{
    // We don't know the ConstraintIndex yet since this hasn't been adopted 
    // by the MatterSubsystem.
    updImpl().updCustomImpl().setMyMatterSubsystem(matter, ConstraintIndex());
}

Constraint::Custom::Implementation::Implementation
   (SimbodyMatterSubsystem& matter, int mp, int mv, int ma) 
:   PIMPLHandle<Implementation,ImplementationImpl>
        (new ImplementationImpl(new CustomImpl(mp,mv,ma))) 
{
    // We don't know the ConstraintIndex yet since this hasn't been adopted 
    // by the MatterSubsystem.
    updImpl().updCustomImpl().setMyMatterSubsystem(matter, ConstraintIndex());
}

const SimbodyMatterSubsystem& 
Constraint::Custom::Implementation::getMatterSubsystem() const {
    return getImpl().getCustomImpl().getMyMatterSubsystem();
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
getOneQ(const State&                                state,
             const Array_<Real,ConstrainedQIndex>&  constrainedQ,
             ConstrainedMobilizerIndex              mobilizer, 
             MobilizerQIndex                        whichQ) const 
{   return getImpl().getCustomImpl()
            .getOneQ(state,constrainedQ,mobilizer,whichQ); }

Real Constraint::Custom::Implementation::
getOneQDot(const State&                                state,
                const Array_<Real,ConstrainedQIndex>&  constrainedQDot,
                ConstrainedMobilizerIndex              mobilizer, 
                MobilizerQIndex                        whichQ) const
{   return getImpl().getCustomImpl()
            .getOneQDot(state,constrainedQDot,mobilizer,whichQ); }

Real Constraint::Custom::Implementation::
getOneQDotDot(const State&                                state,
                   const Array_<Real,ConstrainedQIndex>&  constrainedQDotDot,
                   ConstrainedMobilizerIndex              mobilizer, 
                   MobilizerQIndex                        whichQ) const
{   return getImpl().getCustomImpl()
            .getOneQDotDot(state,constrainedQDotDot,mobilizer,whichQ); }


Real Constraint::Custom::Implementation::
getOneU(const State&                                state,
             const Array_<Real,ConstrainedUIndex>&  constrainedU,
             ConstrainedMobilizerIndex              mobilizer, 
             MobilizerUIndex                        whichU) const
{   return getImpl().getCustomImpl()
            .getOneU(state,constrainedU,mobilizer,whichU); }

Real Constraint::Custom::Implementation::
getOneUDot(const State&                                state,
                const Array_<Real,ConstrainedUIndex>&  constrainedUDot,
                ConstrainedMobilizerIndex              mobilizer, 
                MobilizerUIndex                        whichU) const
{   return getImpl().getCustomImpl()
            .getOneUDot(state,constrainedUDot,mobilizer,whichU); }

Real Constraint::Custom::Implementation::
getOneQFromState(const State& s, ConstrainedMobilizerIndex cmx, MobilizerQIndex mqx) const {
    return getImpl().getCustomImpl().getOneQFromState(s,cmx,mqx);
}

Real Constraint::Custom::Implementation::
getOneUFromState(const State& s, ConstrainedMobilizerIndex cmx, MobilizerUIndex mux) const {
    return getImpl().getCustomImpl().getOneUFromState(s,cmx,mux);
}


Real Constraint::Custom::Implementation::
getOneQDotFromState(const State&                s, 
                    ConstrainedMobilizerIndex   cmx, 
                    MobilizerQIndex             mqx) const {
    return getImpl().getCustomImpl().getOneQDotFromState(s,cmx,mqx);
}


// Apply a generalized (mobility) force to a particular mobility of the given constrained body B,
// adding it in to the appropriate slot of the mobilityForces vector.
void Constraint::Custom::Implementation::
addInOneMobilityForce
   (const State& s, ConstrainedMobilizerIndex M, MobilizerUIndex whichU,
    Real fu, Array_<Real,ConstrainedUIndex>& mobilityForces) const 
{
    getImpl().getCustomImpl().addInOneMobilityForce(s,M,whichU,fu,mobilityForces);
}

void Constraint::Custom::Implementation::
addInOneQForce
   (const State&                    state, 
    ConstrainedMobilizerIndex       mobilizer, 
    MobilizerQIndex                 whichQ,
    Real                            fq, 
    Array_<Real,ConstrainedQIndex>& qForces) const
{
    getImpl().getCustomImpl().addInOneQForce(state,mobilizer,whichQ,fq,qForces);
}


const Transform& Constraint::Custom::Implementation::
getBodyTransform
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getImpl().getCustomImpl().getBodyTransform(allX_AB,bodyB); }

const SpatialVec& Constraint::Custom::Implementation::
getBodyVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getImpl().getCustomImpl().getBodyVelocity(allV_AB,bodyB); }

const SpatialVec& Constraint::Custom::Implementation::
getBodyAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getImpl().getCustomImpl().getBodyAcceleration(allA_AB,bodyB); }

const Transform&  Constraint::Custom::Implementation::
getBodyTransformFromState(const State& s, ConstrainedBodyIndex B) const
{
    return getImpl().getCustomImpl().getBodyTransformFromState(s,B);
}

const SpatialVec& Constraint::Custom::Implementation::
getBodyVelocityFromState(const State& s, ConstrainedBodyIndex B) const
{
    return getImpl().getCustomImpl().getBodyVelocityFromState(s,B);
}

void Constraint::Custom::Implementation::
addInStationForce
   (const State& s, ConstrainedBodyIndex B, const Vec3& p_B, 
    const Vec3& forceInA, 
    Array_<SpatialVec,ConstrainedBodyIndex>& bodyForcesInA) const
{
    getImpl().getCustomImpl().addInStationForce(s,B,p_B,forceInA,bodyForcesInA);
}

void Constraint::Custom::Implementation::
addInBodyTorque
   (const State& s, ConstrainedBodyIndex B,
    const Vec3& torqueInA, 
    Array_<SpatialVec,ConstrainedBodyIndex>& bodyForcesInA) const
{
    getImpl().getCustomImpl().addInBodyTorque(s,B,torqueInA,bodyForcesInA);
}

void Constraint::Custom::Implementation::
getMultipliers(const State&  s, 
               Array_<Real>& multipliers) const
{
    int mp, mv, ma;
    const Constraint::CustomImpl& cimpl = getImpl().getCustomImpl();
    cimpl.getNumConstraintEquationsInUse(s,mp,mv,ma);
    multipliers.resize(mp+mv+ma);
    cimpl.getMultipliers(s, multipliers.size(), &multipliers[0]);
}



// Default implementations for ConstraintImpl virtuals throw "unimplemented"
// exceptions. These shouldn't be called unless the concrete constraint has
// given a non-zero value for mp, mv, and/or ma which is a promise to 
// implement the associated methods.

    // These must be defined if there are any position (holonomic) constraints 
    // defined.

void Constraint::Custom::Implementation::
calcPositionErrors     
   (const State&                                    state,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "calcPositionErrors");
}

void Constraint::Custom::Implementation::
calcPositionDotErrors      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "realizePositionDotErrors");
}

void Constraint::Custom::Implementation::
calcPositionDotDotErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "calcPositionDotDotErrors");
}

void Constraint::Custom::Implementation::
addInPositionConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,ConstrainedQIndex>&             qForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "addInPositionConstraintForces");
}

    // These must be defined if there are any velocity (nonholonomic) 
    // constraints defined.

void Constraint::Custom::Implementation::
calcVelocityErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "calcVelocityErrors");
}


void Constraint::Custom::Implementation::
calcVelocityDotErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "calcVelocityDotErrors");
}


void Constraint::Custom::Implementation::
addInVelocityConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,ConstrainedUIndex>&             mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "addInVelocityConstraintForces");
}



    // These must be defined if there are any acceleration-only constraints 
    // defined.

void Constraint::Custom::Implementation::
calcAccelerationErrors      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::Custom::Implementation", "calcAccelerationErrors");
}

void Constraint::Custom::Implementation::
addInAccelerationConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,ConstrainedUIndex>&             mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
      "Constraint::Custom::Implementation", "addInAccelerationConstraintForces");
}



//==============================================================================
//                       CONSTRAINT::COORDINATE COUPLER
//==============================================================================
Constraint::CoordinateCoupler::CoordinateCoupler
   (SimbodyMatterSubsystem&             matter, 
    const Function*                     function, 
    const Array_<MobilizedBodyIndex>&   coordBody, 
    const Array_<MobilizerQIndex>&      coordIndex)
:   Custom(new CoordinateCouplerImpl(matter, function, coordBody, coordIndex)) 
{}

Constraint::CoordinateCouplerImpl::CoordinateCouplerImpl
   (SimbodyMatterSubsystem&             matter, 
    const Function*                     function, 
    const Array_<MobilizedBodyIndex>&   coordMobod, 
    const Array_<MobilizerQIndex>&      coordQIndex)
:   Implementation(matter, 1, 0, 0), function(function), 
    coordBodies(coordMobod.size()), coordIndices(coordQIndex),
    temp(coordBodies.size()), referenceCount(new int[1]) 
{
    assert(coordBodies.size() == coordIndices.size());
    assert(coordIndices.size() == function->getArgumentSize());
    assert(function->getMaxDerivativeOrder() >= 2);
    referenceCount[0] = 1;
    for (int i = 0; i < (int)coordBodies.size(); ++i) {
        const MobilizedBody& mobod = matter.getMobilizedBody(coordMobod[i]);
        coordBodies[i] =  addConstrainedMobilizer(mobod);
    }
}

void Constraint::CoordinateCouplerImpl::
calcPositionErrors     
   (const State&                                    s,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr) const
{
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQ(s, constrainedQ, coordBodies[i], coordIndices[i]);
    perr[0] = function->calcValue(temp);
}

void Constraint::CoordinateCouplerImpl::
calcPositionDotErrors      
   (const State&                                    s,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) const
{
    pverr[0] = 0;
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQFromState(s, coordBodies[i], coordIndices[i]);
    Array_<int> components(1);
    for (int i = 0; i < temp.size(); ++i) {
        components[0] = i;
        pverr[0] += function->calcDerivative(components, temp)
                    * getOneQDot(s, constrainedQDot, 
                                 coordBodies[i], coordIndices[i]);
    }
}

void Constraint::CoordinateCouplerImpl::
calcPositionDotDotErrors     
   (const State&                                    s,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) const
{
    paerr[0] = 0.0;
    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQFromState(s, coordBodies[i], coordIndices[i]);

    // TODO this could be made faster by using symmetry.
    Array_<int> components(2);
    for (int i = 0; i < temp.size(); ++i) {
        components[0] = i;
        Real qdoti = getOneQDotFromState(s, coordBodies[i], coordIndices[i]);
        for (int j = 0; j < temp.size(); ++j) {
            components[1] = j;
            Real qdotj = getOneQDotFromState(s, coordBodies[j], coordIndices[j]);
            paerr[0] += function->calcDerivative(components, temp)
                        * qdoti * qdotj;
        }
    }

    Array_<int> component(1);
    for (int i = 0; i < temp.size(); ++i) {
        component[0] = i;
        paerr[0] += function->calcDerivative(component, temp)
                    * getOneQDotDot(s, constrainedQDotDot, 
                                    coordBodies[i], coordIndices[i]);
    }
}

void Constraint::CoordinateCouplerImpl::
addInPositionConstraintForces
   (const State&                                s, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForces,
    Array_<Real,ConstrainedQIndex>&             qForces) const
{
    assert(multipliers.size() == 1);
    assert(bodyForces.size() == 0);

    const Real lambda = multipliers[0];

    for (int i = 0; i < temp.size(); ++i)
        temp[i] = getOneQFromState(s, coordBodies[i], coordIndices[i]);

    Array_<int> components(1);
    for (int i = 0; i < temp.size(); ++i) {
        components[0] = i;
        const Real fq = lambda * function->calcDerivative(components, temp);
        addInOneQForce(s, coordBodies[i], coordIndices[i], fq, qForces);
    }
}



//==============================================================================
//                          CONSTRAINT::SPEED COUPLER
//==============================================================================
Constraint::SpeedCoupler::SpeedCoupler
   (SimbodyMatterSubsystem&             matter, 
    const Function*                     function, 
    const Array_<MobilizedBodyIndex>&   speedBody, 
    const Array_<MobilizerUIndex>&      speedIndex)
:   Custom(new SpeedCouplerImpl(matter, function, speedBody, speedIndex, 
                                Array_<MobilizedBodyIndex>(), 
                                Array_<MobilizerQIndex>())) {}

Constraint::SpeedCoupler::SpeedCoupler
   (SimbodyMatterSubsystem&             matter, 
    const Function*                     function, 
    const Array_<MobilizedBodyIndex>&   speedBody, 
    const Array_<MobilizerUIndex>&      speedIndex,
    const Array_<MobilizedBodyIndex>&   coordBody, 
    const Array_<MobilizerQIndex>&      coordIndex)
:   Custom(new SpeedCouplerImpl(matter, function, speedBody, speedIndex, 
                                coordBody, coordIndex)) {}

Constraint::SpeedCouplerImpl::SpeedCouplerImpl
   (SimbodyMatterSubsystem& matter, 
    const Function* function, 
    const Array_<MobilizedBodyIndex>& speedBody, 
    const Array_<MobilizerUIndex>& speedIndex,
    const Array_<MobilizedBodyIndex>& coordBody, 
    const Array_<MobilizerQIndex>& coordIndex)
:   Implementation(matter, 0, 1, 0), function(function), 
    speedBodies(speedBody.size()), speedIndices(speedIndex), 
    coordBodies(coordBody), coordIndices(coordIndex),
    temp(speedBody.size()+coordBody.size()), referenceCount(new int[1]) 
{
    assert(speedBodies.size() == speedIndices.size());
    assert(coordBodies.size() == coordIndices.size());
    assert(temp.size() == function->getArgumentSize());
    assert(function->getMaxDerivativeOrder() >= 2);

    referenceCount[0] = 1;
    for (int i = 0; i < (int)speedBodies.size(); ++i) {
        const MobilizedBody& mobod = matter.getMobilizedBody(speedBody[i]);
        speedBodies[i] = addConstrainedMobilizer(mobod);
    }
}

// Constraint is f(q,u)=0, i.e. verr=f(q,u).
void Constraint::SpeedCouplerImpl::
calcVelocityErrors     
   (const State&                                    s,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr) const
{
    for (int i = 0; i < (int) speedBodies.size(); ++i)
        temp[i] = getOneU(s, constrainedU, speedBodies[i], speedIndices[i]);
    for (int i = 0; i < (int) coordBodies.size(); ++i)
        temp[i+speedBodies.size()] = 
            getMatterSubsystem().getMobilizedBody(coordBodies[i])
                                .getOneQ(s, coordIndices[i]);
    verr[0] = function->calcValue(temp);
}

// d verr / dt = (df/du)*udot + (df/dq)*qdot.
void Constraint::SpeedCouplerImpl::
calcVelocityDotErrors     
   (const State&                                    s,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr) const 
{
    for (int i = 0; i < (int)speedBodies.size(); ++i)
        temp[i] = getOneUFromState(s, speedBodies[i], speedIndices[i]);
    for (int i = 0; i < (int)coordBodies.size(); ++i) {
        const Real q = getMatterSubsystem().getMobilizedBody(coordBodies[i])
                                           .getOneQ(s, coordIndices[i]);
        temp[i+speedBodies.size()] = q;
    }

    Array_<int> components(1);
    vaerr[0] = 0;
    // Differentiate the u-dependent terms here.
    for (int i = 0; i < (int)speedBodies.size(); ++i) {
        components[0] = i;
        vaerr[0] += function->calcDerivative(components, temp)
                    * getOneUDot(s, constrainedUDot, 
                                 speedBodies[i], speedIndices[i]);
    }
    // Differentiate the q-dependent terms here.
    for (int i = 0; i < (int)coordBodies.size(); ++i) {
        components[0] = i + speedBodies.size();
        const Real qdot = getMatterSubsystem().getMobilizedBody(coordBodies[i])
                                              .getOneQDot(s, coordIndices[i]);
        vaerr[0] += function->calcDerivative(components, temp) * qdot;
    }
}

// Force is (df/du)*lambda.
void Constraint::SpeedCouplerImpl::
addInVelocityConstraintForces
   (const State&                                s, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForces,
    Array_<Real,ConstrainedUIndex>&             mobilityForces) const
{
    assert(multipliers.size() == 1);
    const Real lambda = multipliers[0];

    for (int i = 0; i < (int) speedBodies.size(); ++i)
        temp[i] = getOneUFromState(s, speedBodies[i], speedIndices[i]);
    for (int i = 0; i < (int) coordBodies.size(); ++i)
        temp[i+speedBodies.size()] = 
            getMatterSubsystem().getMobilizedBody(coordBodies[i])
                                .getOneQ(s, coordIndices[i]);

    Array_<int> components(1);
    // Only the u-dependent terms generate forces.
    for (int i = 0; i < (int) speedBodies.size(); ++i) {
        components[0] = i;
        const Real force = function->calcDerivative(components, temp)
                           * lambda;
        addInOneMobilityForce(s, speedBodies[i], speedIndices[i], 
                              force, mobilityForces);
    }
}



//==============================================================================
//                        CONSTRAINT::PRESCRIBED MOTION
//==============================================================================
Constraint::PrescribedMotion::PrescribedMotion
   (SimbodyMatterSubsystem&     matter, 
    const Function*             function, 
    MobilizedBodyIndex          coordBody, 
    MobilizerQIndex             coordIndex)
:   Custom(new PrescribedMotionImpl(matter, function, coordBody, coordIndex)) {}

Constraint::PrescribedMotionImpl::PrescribedMotionImpl
   (SimbodyMatterSubsystem& matter, 
    const Function* function, 
    MobilizedBodyIndex coordBody, 
    MobilizerQIndex coordIndex)
:   Implementation(matter, 1, 0, 0), function(function), 
    coordIndex(coordIndex), temp(1), referenceCount(new int[1]) 
{
    assert(function->getArgumentSize() == 1);
    assert(function->getMaxDerivativeOrder() >= 2);

    referenceCount[0] = 1;
    const MobilizedBody& mobod = matter.getMobilizedBody(coordBody);
    this->coordBody = addConstrainedMobilizer(mobod);
}

void Constraint::PrescribedMotionImpl::
calcPositionErrors     
   (const State&                                    s,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr) const
{
    temp[0] = s.getTime();
    perr[0] = getOneQ(s, constrainedQ, coordBody, coordIndex) 
              - function->calcValue(temp);
}

void Constraint::PrescribedMotionImpl::
calcPositionDotErrors      
   (const State&                                    s,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) const
{
    temp[0] = s.getTime();
    Array_<int> components(1, 0); // i.e., components={0}
    pverr[0] = getOneQDot(s, constrainedQDot, coordBody, coordIndex) 
               - function->calcDerivative(components, temp);
}

void Constraint::PrescribedMotionImpl::
calcPositionDotDotErrors     
   (const State&                                    s,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) const
{
    temp[0] = s.getTime();
    Array_<int> components(2, 0); // i.e., components={0,0}
    paerr[0] = getOneQDotDot(s, constrainedQDotDot, coordBody, coordIndex)  
               - function->calcDerivative(components, temp);
}

void Constraint::PrescribedMotionImpl::
addInPositionConstraintForces
   (const State&                                s, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForces,
    Array_<Real,ConstrainedQIndex>&             qForces) const
{
    const Real fq = multipliers[0];
    addInOneQForce(s, coordBody, coordIndex, fq, qForces);
}



//==============================================================================
//                              CONSTRAINT IMPL
//==============================================================================


// =============================================================================
//                              REALIZE TOPOLOGY
// =============================================================================
void ConstraintImpl::realizeTopology(State& s) const {
    // Calculate the relevant Subtree. There might not be any Constrained 
    // Bodies here but we want to make sure we have a properly initialized 
    // empty Subtree in that case.
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



// =============================================================================
//                               REALIZE MODEL
// =============================================================================
void ConstraintImpl::realizeModel(State& s) const
{
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "ConstraintImpl::realizeModel() can't be called until after realizeToplogy().");

    realizeModelVirtual(s); // delegate to concrete constraint
}



// =============================================================================
//                              REALIZE INSTANCE
// =============================================================================
// There are several tasks here that can be performed now that we have values
// for the Instance stage state variables (including the enable/disable flag
// for Constraints):
// (1) Count up the number of holonomic, nonholonomic, and acceleration-only 
//     constraint equations to be contributed by each Constraint, and assign 
//     corresponding slots in constraint-equation ordered arrays, such as the 
//     State's constraint error arrays.
// (2) Count up the number of constrained bodies, mobilizers, and corresoponding
//     constrained mobilities to be affected by each Constraint, and assign
//     corresponding slots for use in "pools" that are indexed by these.
// (3) Above we assigned q's and u's to each mobilizer and stored the results in
//     the Model cache, now we can determine which of those q's and u's are 
//     involved in each constraint. We need to collect up both the set of 
//     directly-constrained q's and u's resulting from ConstrainedMobilizers, 
//     and indirectly-constrained ones arising from their effects on 
//     ConstrainedBodies. Together we call those "participating q's" and 
//     "participating u's" (or "participating mobilities").
// The results of these computations goes in the Instance cache.

void ConstraintImpl::realizeInstance(const State& s) const {
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    const SBInstanceVars& instanceVars  = matter.getInstanceVars(s);
    SBInstanceCache&      ic = matter.updInstanceCache(s);

    const int ncb = getNumConstrainedBodies();
    const int ncm = getNumConstrainedMobilizers();

    // Note that there is a "per constraint info" object for every declared
    // constraint, whether enabled or not.
    SBInstancePerConstraintInfo& cInfo =
        ic.updConstraintInstanceInfo(myConstraintIndex);

    cInfo.clear();
    cInfo.allocateConstrainedMobilizerInstanceInfo(ncm);

    // We're in the process of counting up constraint equations in the
    // totalN...Constraints variables; on entry they are set to the number
    // we've seen so far and we'll increment them here to add in the 
    // contributions from this Constraint.

    if (isDisabled(s)) {
        // Just to be neat, we'll assign zero-width segments where our slots
        // would have gone if this constraint were enabled.
        cInfo.holoErrSegment    = 
            Segment(0,ic.totalNHolonomicConstraintEquationsInUse);
        cInfo.nonholoErrSegment = 
            Segment(0,ic.totalNNonholonomicConstraintEquationsInUse);
        cInfo.accOnlyErrSegment = 
            Segment(0,ic.totalNAccelerationOnlyConstraintEquationsInUse);
        
        cInfo.consBodySegment      = Segment(0, ic.totalNConstrainedBodiesInUse);
        cInfo.consMobilizerSegment = Segment(0, ic.totalNConstrainedMobilizersInUse);
        cInfo.consQSegment = Segment(0, ic.totalNConstrainedQInUse);
        cInfo.consUSegment = Segment(0, ic.totalNConstrainedUInUse);    
        return;
    }

    // This constraint is enabled.

    // These count just the primary contraint equations, not their time 
    // derivatives.
    int mHolo, mNonholo, mAccOnly;
    calcNumConstraintEquationsInUse(s, mHolo, mNonholo, mAccOnly);

    // Must allocate space for primary constraint equations & time derivatives.
    //                                length         offset
    cInfo.holoErrSegment    = Segment(mHolo,    ic.totalNHolonomicConstraintEquationsInUse);
    cInfo.nonholoErrSegment = Segment(mNonholo, ic.totalNNonholonomicConstraintEquationsInUse);
    cInfo.accOnlyErrSegment = Segment(mAccOnly, ic.totalNAccelerationOnlyConstraintEquationsInUse);

    ic.totalNHolonomicConstraintEquationsInUse        += mHolo;
    ic.totalNNonholonomicConstraintEquationsInUse     += mNonholo;
    ic.totalNAccelerationOnlyConstraintEquationsInUse += mAccOnly;  

    cInfo.consBodySegment      = Segment(ncb, ic.totalNConstrainedBodiesInUse);
    cInfo.consMobilizerSegment = Segment(ncm, ic.totalNConstrainedMobilizersInUse);
    ic.totalNConstrainedBodiesInUse     += ncb;
    ic.totalNConstrainedMobilizersInUse += ncm;

    // At this point we can find out how many q's and u's are associated with
    // each of the constrained mobilizers. We'll create packed arrays of q's and
    // u's ordered corresponding to the ConstrainedMobilizerIndices. We'll 
    // record these in the InstanceCache, by storing the ConstrainedQIndex and 
    // ConstrainedUIndex of the lowest-numbered coordinate and mobility 
    // associated with each of the ConstrainedMobilizers, along with the number 
    // of q's and u's.

    for (ConstrainedMobilizerIndex cmx(0); cmx < ncm; ++cmx) {
        SBInstancePerConstrainedMobilizerInfo& mInfo = 
            cInfo.updConstrainedMobilizerInstanceInfo(cmx);

        const MobilizedBodyIndex mbx = 
            getMobilizedBodyIndexOfConstrainedMobilizer(cmx);
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

    // Now we can assign slots for Qs and Us.
    const int ncq = cInfo.getNumConstrainedQ();
    const int ncu = cInfo.getNumConstrainedU();
    cInfo.consQSegment = Segment(ncq, ic.totalNConstrainedQInUse);
    cInfo.consUSegment = Segment(ncu, ic.totalNConstrainedUInUse);
    ic.totalNConstrainedQInUse += ncq;
    ic.totalNConstrainedUInUse += ncu;

    // Now collect all the participating mobilities. This includes the 
    // constrained mobilities as well as every q and u that can affect the 
    // constraint equations which involve constrained bodies. At the end we'll 
    // sort this list by subsystem QIndex/UIndex and remove duplicates.
    cInfo.participatingQ = cInfo.constrainedQ;
    cInfo.participatingU = cInfo.constrainedU;

    const Array_<MobilizedBodyIndex>& bodies = mySubtree.getAllBodies();
    for (int b=1; b<(int)bodies.size(); ++b) { // skip the Ancestor body 0
        QIndex qix; int nq;
        UIndex uix; int nu;
        matter.findMobilizerQs(s,bodies[b],qix,nq);
        matter.findMobilizerUs(s,bodies[b],uix,nu);
        for (int i=0; i<nq; ++i) cInfo.participatingQ.push_back(QIndex(qix+i));
        for (int i=0; i<nu; ++i) cInfo.participatingU.push_back(UIndex(uix+i));
    }

    // Caution: std::unique does not automatically shorten the original list.
    std::sort(cInfo.participatingQ.begin(), cInfo.participatingQ.end());
    Array_<QIndex>::iterator newEnd =
        std::unique(cInfo.participatingQ.begin(), cInfo.participatingQ.end());
    cInfo.participatingQ.erase(newEnd, cInfo.participatingQ.end());

    realizeInstanceVirtual(s); // delegate to concrete constraint
}



// =============================================================================
//                                REALIZE TIME
// =============================================================================
void ConstraintImpl::realizeTime(const SBStateDigest& sbs) const {
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return;

    realizeTimeVirtual(sbs.getState()); // nothing to do in the base class
}



// =============================================================================
//                             REALIZE POSITION
// =============================================================================
void ConstraintImpl::realizePosition(const SBStateDigest& sbs) const {
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return;
    realizePositionVirtual(sbs.getState()); // delegate to concrete constraint
}



// =============================================================================
//                              REALIZE VELOCITY
// =============================================================================
void ConstraintImpl::realizeVelocity(const SBStateDigest& sbs) const {
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return;
    realizeVelocityVirtual(sbs.getState()); // delegate to concrete constraint
}



// =============================================================================
//                              REALIZE DYNAMICS
// =============================================================================
void ConstraintImpl::realizeDynamics(const SBStateDigest& sbs) const {
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return;
    realizeDynamicsVirtual(sbs.getState()); // delegate to concrete constraint
}



// =============================================================================
//                            REALIZE ACCELERATION
// =============================================================================
void ConstraintImpl::realizeAcceleration(const SBStateDigest& sbs) const {
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return; 
    realizeAccelerationVirtual(sbs.getState()); // delegate to concrete constraint
}



// =============================================================================
//                             REALIZE REPORT
// =============================================================================
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

Real ConstraintImpl::getOneQFromState
   (const State& s, ConstrainedMobilizerIndex cmx, MobilizerQIndex whichQ) const
{
    const QIndex qx = getQIndexOfConstrainedQ(s, getConstrainedQIndex(s, cmx, whichQ));
    return getMyMatterSubsystemRep().getQ(s)[qx];
}

Real ConstraintImpl::getOneUFromState
   (const State& s, ConstrainedMobilizerIndex cmx, MobilizerUIndex whichU) const 
{
    const UIndex ux = getUIndexOfConstrainedU(s, getConstrainedUIndex(s, cmx, whichU));
    return getMyMatterSubsystemRep().getU(s)[ux];
}

Real ConstraintImpl::getOneQDotFromState
   (const State&                s, 
    ConstrainedMobilizerIndex   cmx, 
    MobilizerQIndex             whichQ) const
{
    const QIndex qx = getQIndexOfConstrainedQ
                                (s, getConstrainedQIndex(s, cmx, whichQ));
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return matter.getQDot(s)[qx];
}

Real ConstraintImpl::getOneUDotFromState
   (const State&                s,
    ConstrainedMobilizerIndex   cmx, 
    MobilizerUIndex             whichU) const
{
    const UIndex ux = getUIndexOfConstrainedU
                                (s, getConstrainedUIndex(s, cmx, whichU));
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return matter.getUDot(s)[ux];
}


Real ConstraintImpl::getOneQDotDotFromState
   (const State&                s, 
    ConstrainedMobilizerIndex   cmx, 
    MobilizerQIndex             whichQ) const
{
    const QIndex qx = getQIndexOfConstrainedQ(s, getConstrainedQIndex(s, cmx, whichQ));
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    return matter.getQDotDot(s)[qx];
}

// These are used to retrieve the indicated values from the State cache.
const Transform& ConstraintImpl::getBodyTransformFromState
   (const State& s, ConstrainedBodyIndex B) const 
{
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];

    if (!myAncestorBodyIsNotGround) 
        return matter.getBodyTransform(s,bodyB); // X_AB==X_GB

    static const Transform X_AA; // identity Transform
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();
    if (bodyB == ancestorA)
        return X_AA;

    const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[B];
    return matter.getTreePositionCache(s)
                    .constrainedBodyConfigInAncestor[bx]; // X_AB
}
const SpatialVec& ConstraintImpl::getBodyVelocityFromState
   (const State& s, ConstrainedBodyIndex B) const
{
    const SimbodyMatterSubsystemRep& matter    = getMyMatterSubsystemRep();
    const MobilizedBodyIndex         bodyB     = myConstrainedBodies[B];

    if (!myAncestorBodyIsNotGround) 
        return matter.getBodyVelocity(s,bodyB); // V_AB==V_GB

    static const SpatialVec V_AA(Vec3(0),Vec3(0)); // zero velocity
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();
    if (bodyB == ancestorA)
        return V_AA;

    const AncestorConstrainedBodyPoolIndex bx = myPoolIndex[B];
    return matter.getTreeVelocityCache(s)
                    .constrainedBodyVelocityInAncestor[bx]; // V_AB
}



// =============================================================================
//                 CALC CONSTRAINED BODY TRANSFORM IN ANCESTOR
// =============================================================================
// 63 flops per constrained body
void ConstraintImpl::calcConstrainedBodyTransformInAncestor      // X_AB
   (const SBStateDigest& sbs, SBTreePositionCache& tpc) const 
{
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return;
    if (!myAncestorBodyIsNotGround) return;
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();

    // We expect Ground-relative kinematics already to have been calculated in 
    // tpc, but we can't verify that.
    const Transform& X_GA = tpc.getX_GB(ancestorA);

    for (ConstrainedBodyIndex cbx(0); cbx < getNumConstrainedBodies(); ++cbx) {
        const MobilizedBodyIndex bodyB = myConstrainedBodies[cbx];
        if (bodyB == ancestorA) continue; // skip ancestor if it's a constrained body also
        const AncestorConstrainedBodyPoolIndex px = myPoolIndex[cbx];
        assert(px.isValid());

        const Transform& X_GB = tpc.getX_GB(bodyB);
        tpc.updX_AB(px) = ~X_GA*X_GB;  // 63 flops
    }
}



// =============================================================================
//                  CALC CONSTRAINED BODY VELOCITY IN ANCESTOR
// =============================================================================
// 51 flops per constrained body
void ConstraintImpl::calcConstrainedBodyVelocityInAncestor       // V_AB
   (const SBStateDigest& sbs, SBTreeVelocityCache& tvc) const 
{
    const SBInstanceVars& instanceVars  = sbs.getInstanceVars();
    if (instanceVars.disabled[myConstraintIndex]) return;
    if (!myAncestorBodyIsNotGround) return;
    const MobilizedBodyIndex ancestorA = mySubtree.getAncestorMobilizedBodyIndex();

    const SBTreePositionCache& tpc = sbs.getTreePositionCache();

    // All position kinematics has been calculated, and we also expect 
    // Ground-relative velocity kinematics already to have been calculated in tvc, 
    // but we can't verify that.
    const Transform&  X_GA = tpc.getX_GB(ancestorA);
    const SpatialVec& V_GA = tvc.getV_GB(ancestorA); 

    for (ConstrainedBodyIndex cbx(0); cbx < getNumConstrainedBodies(); ++cbx) {
        const MobilizedBodyIndex bodyB = myConstrainedBodies[cbx];
        if (bodyB == ancestorA) continue; // skip the ancestor itself if it is a constrained body also
        const AncestorConstrainedBodyPoolIndex px = myPoolIndex[cbx];
        assert(px.isValid());

        const Transform&  X_GB = tpc.getX_GB(bodyB);
        const SpatialVec& V_GB = tvc.getV_GB(bodyB);

        // 6 flops
        const Vec3 p_AB_G     = X_GB.p() - X_GA.p();
        const Vec3 p_AB_G_dot = V_GB[1]  - V_GA[1];        // d/dt p taken in G

        // 3 flops
        const Vec3 w_AB_G = V_GB[0] - V_GA[0];             // relative angular velocity of B in A, exp. in G

        // To get d/dt p taken in A, get derivative in G and remove the contribution generated by
        // A's velocity in G.
        // 12 flops
        const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // time deriv of p in A, exp in G

        // 30 flops
        tvc.updV_AB(px) = ~X_GA.R() * SpatialVec(w_AB_G, v_AB_G);     // re-express in A
    }
}


// Find out how many holonomic (position), nonholonomic (velocity),
// and acceleration-only constraint equations are generated by this Constraint
// as it is currently being modeled.
void ConstraintImpl::getNumConstraintEquationsInUse
   (const State& s, int& mHolo, int& mNonholo, int& mAccOnly) const 
{
    const SBInstancePerConstraintInfo& cInfo = 
        getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex);

    mHolo    = cInfo.holoErrSegment.length;
    mNonholo = cInfo.nonholoErrSegment.length;
    mAccOnly = cInfo.accOnlyErrSegment.length;
}

void ConstraintImpl::
getIndexOfMultipliersInUse(const State&     s,
                           MultiplierIndex& px0, 
                           MultiplierIndex& vx0, 
                           MultiplierIndex& ax0) const
{
    const SBInstanceCache& ic = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);

    // Find the offset to our first multiplier in the ModelCache.
    const int firstHoloErr = cInfo.holoErrSegment.offset;
    const int mHolo        = cInfo.holoErrSegment.length;

    const int firstNonholoErr = 
          ic.totalNHolonomicConstraintEquationsInUse //total for whole subsystem
        + cInfo.nonholoErrSegment.offset;
    const int mNonholo = cInfo.nonholoErrSegment.length;

    const int firstAccOnlyErr = 
          ic.totalNHolonomicConstraintEquationsInUse
        + ic.totalNNonholonomicConstraintEquationsInUse
        + cInfo.accOnlyErrSegment.offset;
    const int mAccOnly = cInfo.accOnlyErrSegment.length;

    px0.invalidate(); if (mHolo)    px0=MultiplierIndex(firstHoloErr);
    vx0.invalidate(); if (mNonholo) vx0=MultiplierIndex(firstNonholoErr);
    ax0.invalidate(); if (mAccOnly) ax0=MultiplierIndex(firstAccOnlyErr);
}

void ConstraintImpl::
setMyPartInConstraintSpaceVector(const State& s,
                                 const Vector& myPart,
                                 Vector& constraintSpace) const
{
    const SBInstanceCache& ic = getInstanceCache(s);
    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;

    const int m = mHolo+mNonholo+mAccOnly;
    if (constraintSpace.size() == 0) {
        constraintSpace.resize(m);
        constraintSpace.setToZero();
    }

    SimTK_ERRCHK2_ALWAYS(constraintSpace.size()==m,
        "Constraint::setMyPartInConstraintSpaceVector()",
        "Output vector had size %d but should have had size zero or m=%d.",
        constraintSpace.size(), m);

    int mp, mv, ma;
    getNumConstraintEquationsInUse(s, mp, mv, ma);
    MultiplierIndex px0, vx0, ax0;
    getIndexOfMultipliersInUse(s, px0, vx0, ax0);

    for (int i=0; i<mp; ++i)
        constraintSpace[px0+i] = myPart[i];
    for (int i=0; i<mv; ++i)
        constraintSpace[vx0+i] = myPart[mp+i];
    for (int i=0; i<ma; ++i)
        constraintSpace[ax0+i] = myPart[mp+mv+i];
}

void ConstraintImpl::
getMyPartFromConstraintSpaceVector(const State& s,
                                   const Vector& constraintSpace,
                                   Vector& myPart) const
{
    const SBInstanceCache& ic = getInstanceCache(s);
    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;

    const int m = mHolo+mNonholo+mAccOnly;

    SimTK_ERRCHK2_ALWAYS(constraintSpace.size()==m,
        "Constraint::getMyPartFromConstraintSpaceVector()",
        "Input vector had size %d but should have had size m=%d.",
        constraintSpace.size(), m);

    int mp, mv, ma;
    getNumConstraintEquationsInUse(s, mp, mv, ma);
    MultiplierIndex px0, vx0, ax0;
    getIndexOfMultipliersInUse(s, px0, vx0, ax0);

    myPart.resize(mp+mv+ma);

    for (int i=0; i<mp; ++i)
         myPart[i] = constraintSpace[px0+i];
    for (int i=0; i<mv; ++i)
         myPart[mp+i] = constraintSpace[vx0+i];
    for (int i=0; i<ma; ++i)
         myPart[mp+mv+i] = constraintSpace[ax0+i];
}

void ConstraintImpl::setDisabled(State& s, bool shouldBeDisabled) const {
    getMyMatterSubsystemRep().setConstraintIsDisabled(s, myConstraintIndex, shouldBeDisabled);
}

bool ConstraintImpl::isDisabled(const State& s) const {
    return getMyMatterSubsystemRep().isConstraintDisabled(s, myConstraintIndex);
}

// Call this during construction phase to add a body to the topological 
// structure of this Constraint. This body's mobilizer's mobilities are 
// *not* part of the constraint; mobilizers must be added separately. It is OK
// to add the same body multiple times; it will only get inserted once and 
// you'll get the same index every time.
ConstrainedBodyIndex ConstraintImpl::
addConstrainedBody(const MobilizedBody& b) {
    assert(isInSameSubsystem(b));
    invalidateTopologyCache();

    const ConstrainedBodyIndex nextIx((int)myConstrainedBodies.size());

    // Add to the Mobilized->Constrained map and check for duplicates.
    std::pair<MobilizedBody2ConstrainedBodyMap::iterator, bool> result;
    result = myMobilizedBody2ConstrainedBodyMap.insert(
        MobilizedBody2ConstrainedBodyMap::value_type(b.getMobilizedBodyIndex(), 
                                                     nextIx));
    if (!result.second) {
        // It was already there.
        return result.first->second; // the index we assigned before
    }

    // This is a new constrained body -- add it to the 
    // ConstrainedBody->MobilizedBody map too.
    myConstrainedBodies.push_back(b.getMobilizedBodyIndex());
    return nextIx;
}

// Call this during construction phase to add a mobilizer to the topological 
// structure of this Constraint. All the coordinates q and mobilities u for this
// mobilizer are added also, but we don't know how many of those there will be 
// until Stage::Model. It is OK to add the same mobilizer multiple times; it will
// only get inserted once and you'll get the same index every time. 
ConstrainedMobilizerIndex ConstraintImpl::
addConstrainedMobilizer(const MobilizedBody& b) {
    assert(isInSameSubsystem(b));
    invalidateTopologyCache();

    const ConstrainedMobilizerIndex nextIx((int)myConstrainedMobilizers.size());

    // Add to the Mobilized->Constrained map and check for duplicates.
    std::pair<MobilizedBody2ConstrainedMobilizerMap::iterator, bool> result;
    result = myMobilizedBody2ConstrainedMobilizerMap.insert
       (MobilizedBody2ConstrainedMobilizerMap::value_type
                                            (b.getMobilizedBodyIndex(), nextIx));
    
    if (!result.second) {
        // It was already there.
        return result.first->second; // the index we assigned before
    }
    
    // This is a new constrained mobilizer -- add it to the 
    // ConstrainedMobilizer->MobilizedBody map too.
    myConstrainedMobilizers.push_back(b.getMobilizedBodyIndex());
    return nextIx;
}

QIndex ConstraintImpl::getQIndexOfConstrainedQ(const State& s, ConstrainedQIndex cqx) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const SBInstancePerConstraintInfo& 
        cInfo = ic.getConstraintInstanceInfo(myConstraintIndex);
    return cInfo.getQIndexFromConstrainedQ(cqx);
}

UIndex ConstraintImpl::getUIndexOfConstrainedU(const State& s, ConstrainedUIndex cux) const {
    const SBInstanceCache&             ic    = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);
    return cInfo.getUIndexFromConstrainedU(cux);
}

int ConstraintImpl::getNumConstrainedQ(const State& s) const {
    return getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex)
                              .getNumConstrainedQ();
}

int ConstraintImpl::getNumConstrainedQ
   (const State& s, ConstrainedMobilizerIndex M) const
{
    const SBInstancePerConstrainedMobilizerInfo& mInfo =
        getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex)
                           .getConstrainedMobilizerInstanceInfo(M);
    return mInfo.nQInUse; // same as corresponding Mobod, or 0 if disabled
}

ConstrainedQIndex ConstraintImpl::getConstrainedQIndex
   (const State& s, ConstrainedMobilizerIndex M, MobilizerQIndex which) const 
{
    const int nq = getNumConstrainedQ(s,M);
    assert(0 <= which && which < nq);
    const SBInstancePerConstrainedMobilizerInfo& mInfo =
        getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex)
                           .getConstrainedMobilizerInstanceInfo(M);
    return ConstrainedQIndex(mInfo.firstConstrainedQIndex + which);
}       

int ConstraintImpl::getNumConstrainedU(const State& s) const {
    return getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex)
                              .getNumConstrainedU();
}

int ConstraintImpl::getNumConstrainedU
   (const State& s, ConstrainedMobilizerIndex M) const
{
    const SBInstancePerConstrainedMobilizerInfo& mInfo =
        getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex)
                           .getConstrainedMobilizerInstanceInfo(M);
    return mInfo.nUInUse; // same as corr. MobilizedBody, or 0 if disabled
}

ConstrainedUIndex ConstraintImpl::getConstrainedUIndex
   (const State& s, ConstrainedMobilizerIndex M, MobilizerUIndex which) const 
{
    const int nu = getNumConstrainedU(s,M);
    assert(0 <= which && which < nu);
    const SBInstancePerConstrainedMobilizerInfo& mInfo =
        getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex)
                           .getConstrainedMobilizerInstanceInfo(M);
    return ConstrainedUIndex(mInfo.firstConstrainedUIndex + which);
}  



//==============================================================================
//                       CONVERT Q FORCES TO U FORCES
//==============================================================================
// uForces = ~N * qForces
void ConstraintImpl::
convertQForcesToUForces(const State&                          s, 
                        const Array_<Real,ConstrainedQIndex>& qForces,
                        Array_<Real,ConstrainedUIndex>&       uForces) const
{
    const int ncm = getNumConstrainedMobilizers();
    if (ncm == 0)
        return; // very common!

    const SBInstanceCache& ic = getInstanceCache(s);
    const SBInstancePerConstraintInfo& 
        cInfo = ic.getConstraintInstanceInfo(myConstraintIndex);

    assert(cInfo.getNumConstrainedMobilizers() == ncm);

    const int ncu = cInfo.getNumConstrainedU();
    const int ncq = cInfo.getNumConstrainedQ();

    assert(qForces.size() == ncq);
    uForces.resize(ncu);

    if (ncu == 0) 
        return;

    for (ConstrainedMobilizerIndex cmx(0); cmx < ncm; ++cmx) {
        const SBInstancePerConstrainedMobilizerInfo& mInfo =
            cInfo.getConstrainedMobilizerInstanceInfo(cmx);
        const int nq = mInfo.nQInUse, nu = mInfo.nUInUse;
        if (nu == 0) continue;

        const Real* firstQ = &qForces[mInfo.firstConstrainedQIndex];
        Real* firstU = &uForces[mInfo.firstConstrainedUIndex];
        const ArrayViewConst_<Real,MobilizerQIndex> fq(firstQ, firstQ+nq);
        ArrayView_<Real,MobilizerUIndex> fu(firstU, firstU+nu);
        const MobilizedBody& mobod = 
            getMobilizedBodyFromConstrainedMobilizer(cmx);
        mobod.convertQForceToUForce(s, fq, fu);
    }
}



//==============================================================================
//               CONVERT BODY ACCEL TO CONSTRAINED BODY ACCEL
//==============================================================================
void ConstraintImpl::convertBodyAccelToConstrainedBodyAccel
   (const State&                                    s,
    const Array_<SpatialVec, MobilizedBodyIndex>&   allA_GB,
    Array_<SpatialVec, ConstrainedBodyIndex>&       A_AB) const
{
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    assert(allA_GB.size() == matter.getNumBodies());

    const int ncb = getNumConstrainedBodies();
    A_AB.resize(ncb);

    // If the Ancestor is Ground we're just reordering. 
    if (!isAncestorDifferentFromGround()) {
        for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) {
            const MobilizedBodyIndex mbx = 
                getMobilizedBodyIndexOfConstrainedBody(cbx);
            A_AB[cbx] = allA_GB[mbx];
        }
        return;
    }

    // If the Ancestor is not Ground we have to transform the accelerations 
    // from Ground to Ancestor, at a cost of 105 flops/constrained body (not 
    // just re-expressing).

    const MobilizedBody& ancestor = getAncestorMobilizedBody();
    const MobilizedBodyIndex acx = ancestor.getMobilizedBodyIndex();
    const Transform&  X_GA  = ancestor.getBodyTransform(s);
    const SpatialVec& V_GA  = ancestor.getBodyVelocity(s);
    const SpatialVec& A_GA  = allA_GB[acx]; // from input argument

    for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) {
        const MobilizedBodyIndex mbx =
            getMobilizedBodyIndexOfConstrainedBody(cbx);
        if (mbx == acx) { // common: Ancestor is a constrained body
            A_AB[cbx] = SpatialVec(Vec3(0));
            continue;
        }
        const MobilizedBody& consBody = matter.getMobilizedBody(mbx);
        const Transform&  X_GB = consBody.getBodyTransform(s);
        const SpatialVec& V_GB = consBody.getBodyVelocity(s);
        const SpatialVec& A_GB = allA_GB[mbx];  // from input arg
        A_AB[cbx] = findRelativeAcceleration(X_GA, V_GA, A_GA,
                                             X_GB, V_GB, A_GB);
    }
}



//==============================================================================
//             CONVERT BODY VELOCITY TO CONSTRAINED BODY VELOCITY
//==============================================================================
void ConstraintImpl::convertBodyVelocityToConstrainedBodyVelocity
   (const State&                                    s,
    const Array_<SpatialVec, MobilizedBodyIndex>&   allV_GB,
    Array_<SpatialVec, ConstrainedBodyIndex>&       V_AB) const
{
    const SimbodyMatterSubsystemRep& matter = getMyMatterSubsystemRep();
    assert(allV_GB.size() == matter.getNumBodies());

    const int ncb = getNumConstrainedBodies();
    V_AB.resize(ncb);

    // If the Ancestor is Ground we're just reordering. 
    if (!isAncestorDifferentFromGround()) {
        for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) {
            const MobilizedBodyIndex mbx = 
                getMobilizedBodyIndexOfConstrainedBody(cbx);
            V_AB[cbx] = allV_GB[mbx];
        }
        return;
    }

    // If the Ancestor is not Ground we have to transform the velocities 
    // from Ground to Ancestor, at a cost of 51 flops/constrained body (not 
    // just re-expressing).

    const MobilizedBody& ancestor = getAncestorMobilizedBody();
    const MobilizedBodyIndex acx = ancestor.getMobilizedBodyIndex();
    const Transform&  X_GA  = ancestor.getBodyTransform(s);
    const SpatialVec& V_GA  = allV_GB[acx]; // from input argument

    for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) {
        const MobilizedBodyIndex mbx =
            getMobilizedBodyIndexOfConstrainedBody(cbx);
        if (mbx == acx) { // common: Ancestor is a constrained body
            V_AB[cbx] = SpatialVec(Vec3(0));
            continue;
        }
        const MobilizedBody& consBody = matter.getMobilizedBody(mbx);
        const Transform&  X_GB = consBody.getBodyTransform(s);
        const SpatialVec& V_GB = allV_GB[mbx];  // from input arg
        V_AB[cbx] = findRelativeVelocity(X_GA, V_GA,    // 51 flops
                                         X_GB, V_GB);
    }
}



//==============================================================================
//               GET POSITION, VELOCITY, ACCELERATION ERRORS
//==============================================================================

// Given a state realized to Position stage, extract the position constraint 
// errors corresponding to this Constraint. The 'mp' argument is for sanity 
// checking -- it is an error if that isn't an exact match for the current 
// number of holonomic constraint equations generated by this Constraint. We 
// expect that perr points to an array of at least mp elements that we can 
// write on.
void ConstraintImpl::getPositionErrors(const State& s, int mp, Real* perr) const {
    const SBInstancePerConstraintInfo& cInfo = 
        getInstanceCache(s).getConstraintInstanceInfo(myConstraintIndex);

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
    const SBInstanceCache&             ic    = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);

    assert(mpv ==  cInfo.holoErrSegment.length
                 + cInfo.nonholoErrSegment.length);

    // Get reference to all uerr's for the subsystem.
    const Vector& uerr = getMyMatterSubsystemRep().getUErr(s);

    // Find the offset to our first uerr in the ModelCache.
    const int firstHoloErr = cInfo.holoErrSegment.offset;
    const int mHolo        = cInfo.holoErrSegment.length;

    for (int i=0; i < mHolo; ++i)
        pverr[i] = uerr[firstHoloErr+i];

    const int firstNonholoErr = ic.totalNHolonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.nonholoErrSegment.offset;
    const int mNonholo        = cInfo.nonholoErrSegment.length;

    for (int i=0; i < mNonholo; ++i)
        pverr[mHolo+i] = uerr[firstNonholoErr+i];
}

// Given a State realized to Acceleration stage, extract the acceleration 
// constraint errors corresponding to this Constraint. This includes 
// acceleration constraints which were produced by twice differentiation of 
// holonomic (position) constraints, and differentiation of nonholonomic 
// (velocity) constraints, and acceleration-only constraints which are
// first introduced at the acceleration level. The 'mpva' argument is
// for sanity checking -- it is an error if that isn't an exact match for the
// current number of holonomic+nonholonomic+accelerationOnly (mp+mv+ma) 
// constraint equations generated by this Constraint. We expect that pvaerr 
// points to an array of at least mp+mv+ma elements that we can write on.
void ConstraintImpl::getAccelerationErrors(const State& s, int mpva, Real* pvaerr) const {
    const SBInstanceCache&             ic    = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);

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

    const int firstNonholoErr = ic.totalNHolonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.nonholoErrSegment.offset;
    const int mNonholo        = cInfo.nonholoErrSegment.length;

    for (int i=0; i < mNonholo; ++i)
        pvaerr[mHolo+i] = udoterr[firstNonholoErr+i];

    const int firstAccOnlyErr = ic.totalNHolonomicConstraintEquationsInUse
                                + ic.totalNNonholonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.accOnlyErrSegment.offset;
    const int mAccOnly        = cInfo.accOnlyErrSegment.length;

    for (int i=0; i < mAccOnly; ++i)
        pvaerr[mHolo+mNonholo+i] = udoterr[firstAccOnlyErr+i];
}

// Given a State realized to Acceleration stage, extract the Lagrange
// multipliers corresponding to this Constraint. The 'mpva' argument is for 
// sanity checking -- it is an error if that isn't an exact match for the
// current number of holonomic+nonholonomic+accelerationOnly (mp+mv+ma) 
// constraint equations generated by this Constraint. We expect that lambda 
// points to an array of at least mp+mv+ma elements that we can write on.
void ConstraintImpl::getMultipliers(const State& s, int mpva, Real* lambda) const {
    const SBInstanceCache&             ic    = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);

    assert(mpva ==   cInfo.holoErrSegment.length
                   + cInfo.nonholoErrSegment.length
                   + cInfo.accOnlyErrSegment.length);

    // Get reference to all multipliers for the subsystem. Use "upd" here
    // because we might still be realizing this state.
    const Vector& multipliers = getMyMatterSubsystemRep().updMultipliers(s);

    // Find the offset to our first multiplier in the ModelCache.
    const int firstHoloErr = cInfo.holoErrSegment.offset;
    const int mHolo        = cInfo.holoErrSegment.length;

    for (int i=0; i < mHolo; ++i)
        lambda[i] = multipliers[firstHoloErr+i];

    const int firstNonholoErr = ic.totalNHolonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.nonholoErrSegment.offset;
    const int mNonholo        = cInfo.nonholoErrSegment.length;

    for (int i=0; i < mNonholo; ++i)
        lambda[mHolo+i] = multipliers[firstNonholoErr+i];

    const int firstAccOnlyErr = ic.totalNHolonomicConstraintEquationsInUse
                                + ic.totalNNonholonomicConstraintEquationsInUse // total for whole subsystem
                                + cInfo.accOnlyErrSegment.offset;
    const int mAccOnly        = cInfo.accOnlyErrSegment.length;

    for (int i=0; i < mAccOnly; ++i)
        lambda[mHolo+mNonholo+i] = multipliers[firstAccOnlyErr+i];
}

// Reference this constraint's assigned partition within the larger array.
ArrayView_<SpatialVec,ConstrainedBodyIndex> ConstraintImpl::
updConstrainedBodyForces(const State&        state,
                         Array_<SpatialVec>& allConsBodyForces) const 
{
    const SBInstanceCache& ic = getInstanceCache(state);
    assert(allConsBodyForces.size() == ic.totalNConstrainedBodiesInUse);

    const ConstraintIndex              cx    = getMyConstraintIndex();
    const SBInstancePerConstraintInfo& cInfo = ic.getConstraintInstanceInfo(cx);

    const Segment& consBodySegment = cInfo.consBodySegment;
    const int ncb = consBodySegment.length;

    // No heap allocation is being done here. The longer array can be
    // empty so we're using begin() here rather than &array[0] which 
    // would be illegal in that case. This pointer may be null.
    SpatialVec* firstBodySlot = allConsBodyForces.begin() 
                                + consBodySegment.offset;

    return ArrayView_<SpatialVec,ConstrainedBodyIndex>
                (firstBodySlot, firstBodySlot + ncb);
}

ArrayView_<Real,ConstrainedUIndex> ConstraintImpl::
updConstrainedMobilityForces(const State&  state,
                             Array_<Real>& allConsMobForces) const
{
    const SBInstanceCache& ic = getInstanceCache(state);
    assert(allConsMobForces.size() == ic.totalNConstrainedUInUse);

    const ConstraintIndex              cx    = getMyConstraintIndex();
    const SBInstancePerConstraintInfo& cInfo = ic.getConstraintInstanceInfo(cx);

    const Segment& consUSegment = cInfo.consUSegment;
    const int ncu = consUSegment.length;

    // No heap allocation is being done here. The longer array can be
    // empty so we're using begin() here rather than &array[0] which 
    // would be illegal in that case. This pointer may be null.
    Real* firstMobSlot = allConsMobForces.begin() 
                         + consUSegment.offset;

    return ArrayView_<Real,ConstrainedUIndex>
                (firstMobSlot, firstMobSlot + ncu);
}

const SBInstanceCache& ConstraintImpl::getInstanceCache(const State& s) const {
    return getMyMatterSubsystemRep().getInstanceCache(s);
}
const SBModelCache& ConstraintImpl::getModelCache(const State& s) const {
    return getMyMatterSubsystemRep().getModelCache(s);
}
const SBTreePositionCache& ConstraintImpl::getTreePositionCache(const State& s) const {
    return getMyMatterSubsystemRep().getTreePositionCache(s);
}
const SBTreeVelocityCache& ConstraintImpl::getTreeVelocityCache(const State& s) const {
    return getMyMatterSubsystemRep().getTreeVelocityCache(s);
}
const SBTreeAccelerationCache& ConstraintImpl::getTreeAccelerationCache(const State& s) const {
    return getMyMatterSubsystemRep().getTreeAccelerationCache(s);
}
const SBConstrainedAccelerationCache& ConstraintImpl::
getConstrainedAccelerationCache(const State& s) const {
    return getMyMatterSubsystemRep().getConstrainedAccelerationCache(s);
}
SBConstrainedAccelerationCache& ConstraintImpl::
updConstrainedAccelerationCache(const State& s) const {
    return getMyMatterSubsystemRep().updConstrainedAccelerationCache(s);
}

//------------------------------------------------------------------------------
// Default implementations for ConstraintImpl virtuals throw "unimplemented"
// exceptions. These shouldn't be called unless the concrete constraint has
// given a non-zero value for mp, mv, and/or ma which is a promise to 
// implement the associated methods.
//------------------------------------------------------------------------------

    // These four must be defined if there are any position (holonomic) 
    // constraints defined.

void ConstraintImpl::calcPositionErrorsVirtual      
   (const State&                                    state,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)
    const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "calcPositionErrorsVirtual");
}

void ConstraintImpl::calcPositionDotErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)
    const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "calcPositionDotErrorsVirtual");
}

void ConstraintImpl::calcPositionDotDotErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)
    const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "calcPositionDotDotErrorsVirtual");
}

void ConstraintImpl::addInPositionConstraintForcesVirtual
   (const State&                                    state,
    const Array_<Real>&                             multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForces,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "addInPositionConstraintForcesVirtual");
}


    // These three must be defined if there are any velocity (nonholonomic) 
    // constraints defined.

void ConstraintImpl::calcVelocityErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)
    const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "calcVelocityErrorsVirtual");
}

void ConstraintImpl::calcVelocityDotErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)
    const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "calcVelocityDotErrorsVirtual");
}

void ConstraintImpl::addInVelocityConstraintForcesVirtual
   (const State&                                    state,
    const Array_<Real>&                             multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForces,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "addInVelocityConstraintForcesVirtual");
}


    // These two must be defined if there are any acceleration-only constraints
    // defined.

void ConstraintImpl::calcAccelerationErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr)
    const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "calcAccelerationErrorsVirtual");
}

void ConstraintImpl::addInAccelerationConstraintForcesVirtual
   (const State&                                    state,
    const Array_<Real>&                             multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForces,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "ConstraintImpl", "addInAccelerationConstraintForcesVirtual");
}

//------------------------------------------------------------------------------
// These are interfaces to the constraint operators which first extract
// the operands from a given state. These are thus suitable for use when
// realizing that state at the point where the constraint operator results
// are about to go into the state cache. The cache is not updated here,
// however. Instead the result is returned explicitly in an argument.
//------------------------------------------------------------------------------

// Calculate the mp position errors that would result from the configuration 
// present in the supplied state (that is, q's and body transforms). The state
// must be realized through Time stage and part way through realization of
// Position stage.
void ConstraintImpl::
calcPositionErrorsFromState(const State& s, Array_<Real>& perr) const {
    const SBInstanceCache&             ic    = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);
    const Vector& q = s.getQ();

    const int ncb = getNumConstrainedBodies();
    const int ncq = cInfo.getNumConstrainedQ();

    Array_<Transform, ConstrainedBodyIndex> X_AB(ncb);
    Array_<Real, ConstrainedQIndex>         cq(ncq);

    for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) 
        X_AB[cbx] = getBodyTransformFromState(s, cbx);

    for (ConstrainedQIndex cqx(0); cqx < ncq; ++cqx)
        cq[cqx] = q[cInfo.getQIndexFromConstrainedQ(cqx)];
        
    calcPositionErrors(s,X_AB,cq,perr);
}

// Calculate the mp velocity errors resulting from pdot equations, given a
// configuration and velocities in the supplied state which must be realized
// through Position stage and part way through realization of Velocity stage.
void ConstraintImpl::
calcPositionDotErrorsFromState(const State& s, Array_<Real>& pverr) const {
    const SBInstanceCache&             ic    = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);

    // We're not checking, but the tree velocity cache better have been
    // marked valid by now, indicating that we finished calculating qdots.
    // The state doesn't know about that yet, so we have to use updQDot()
    // here rather than getQDot().
    const Vector& qdot = s.updQDot();

    const int ncb = getNumConstrainedBodies();
    const int ncq = cInfo.getNumConstrainedQ();

    Array_<SpatialVec, ConstrainedBodyIndex> V_AB(ncb);
    Array_<Real, ConstrainedQIndex>          cqdot(ncq);

    for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) 
        V_AB[cbx] = getBodyVelocityFromState(s, cbx);

    for (ConstrainedQIndex cqx(0); cqx < ncq; ++cqx)
        cqdot[cqx] = qdot[cInfo.getQIndexFromConstrainedQ(cqx)];

    calcPositionDotErrors(s,V_AB,cqdot,pverr);
}



// Calculate the mv velocity errors resulting from the nonholonomic constraint
// equations here, taking the configuration and velocities (u, qdot, body
// spatial velocities) from the supplied state, which must be realized through
// Position stage and part way through realization of Velocity stage.
void ConstraintImpl::
calcVelocityErrorsFromState(const State& s, Array_<Real>& verr) const {
    const SBInstanceCache&             ic = getInstanceCache(s);
    const SBInstancePerConstraintInfo& cInfo = 
        ic.getConstraintInstanceInfo(myConstraintIndex);
    const Vector& u = s.getU();

    const int ncb = getNumConstrainedBodies();
    const int ncu = cInfo.getNumConstrainedU();

    Array_<SpatialVec, ConstrainedBodyIndex> V_AB(ncb);
    Array_<Real, ConstrainedUIndex>          cu(ncu);

    for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx) 
        V_AB[cbx] = getBodyVelocityFromState(s, cbx);

    for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
        cu[cux] = u[cInfo.getUIndexFromConstrainedU(cux)];

    calcVelocityErrorsVirtual(s,V_AB,cu,verr);
}



} // namespace SimTK

