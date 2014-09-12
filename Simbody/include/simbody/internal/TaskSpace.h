/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
 * Contributors: Michael Sherman                                              *
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

#include "simbody/internal/common.h"
#include "SimTKmath.h"
#include "SimbodyMatterSubsystem.h"

namespace SimTK {

#define TASKSPACEQUANTITY_GETTER_HELPER(CLASS, MEMVAR, SHORT) \
        realize ## CLASS(s); \
        return Value<CLASS>::downcast( \
        m_matter.getCacheEntry(s, m_ ## MEMVAR ## Index)); \

#define TASKSPACEQUANTITY_MEMBERS(CLASS, MEMVAR, SHORT) \
private: \
    void realize ## CLASS(const State& s) const { \
        if (m_matter.isCacheValueRealized(s, m_ ## MEMVAR ## Index)) \
            return; \
        CLASS& qty = upd ## CLASS(s); \
        if (!qty.m_tspace) qty.m_tspace = this; \
        qty.realize(s); \
        m_matter.markCacheValueRealized(s, m_ ## MEMVAR ## Index); \
    } \
    CLASS& upd ## CLASS(const State& s) const { \
        return Value<CLASS>::updDowncast( \
        m_matter.updCacheEntry(s, m_ ## MEMVAR ## Index)); \
    } \
    CacheEntryIndex m_ ## MEMVAR ## Index; \
public: \
    const CLASS& get ## CLASS(const State& s) const { \
        TASKSPACEQUANTITY_GETTER_HELPER(CLASS, MEMVAR, SHORT); \
    } \
    const CLASS& SHORT(const State& s) const { \
        TASKSPACEQUANTITY_GETTER_HELPER(CLASS, MEMVAR, SHORT); \
    } \

/** (Experimental – API will change – use at your own risk) Efficient and
* convenient computation of quantities necessary for task-space
* model-based controllers.
*
* This class provides convenient access to the quantities used in task-space
* or operational-space controllers, such as the jacobian and the task-space
* mass matrix. Each such quantity is encapsulated in its own class, and objects
* of these types can be used as if they were the quantities (e.g. matrix)
* themselves. This encapsulation allows us to perform
* the necessary calculations more efficiently under the covers.
*
* Computing quantities such as the jacobian \a explicity is often very
* inefficient, and is usually unnecessary. The jacobian usually appears in a
* matrix-vector product, which is more efficient to compute. This class and
* the classes within use operator overloading to allow a user to write code
* as if they had the matrices explicitly; internally, we compute matrix-vector
* products instead. Also, this class caches the quantities once they've
* been computed.
*
* <h3>List of task-space quantities</h3>
*
* In this class, we use Khatib's notation for operational-space controllers [1].
* This notation conflicts with the notation used elsewhere in Simbody.
*  - nu: number of degrees of freedom.
*  - nt: number of tasks (e.g., number of stations).
*  - nst: number of scalar tasks; 3 * nt.
*  - \f$ A \f$ (nu x nu): joint-space mass matrix (M elsewhere in Simbody).
*  - \f$ b \f$ (nu x 1): joint-space inertial forces (Coriolis, etc.).
*  - \f$ g \f$ (nu x 1): joint-space gravity forces.
*  - \f$ \Lambda \f$ (nst x nst): task-space mass matrix.
*  - \f$ p \f$ (nst x 1): task-space gravity forces.
*  - \f$ \mu \f$ (nst x 1): task-space inertial forces.
*  - \f$ J \f$ (nst x nu): task jacobian.
*  - \f$ \bar{J} \f$ (nu x nst): dynamically consistent generalized inverse of the jacobian.
*  - \f$ N \f$ (nu x nu): nullspace projection matrix.
*
*  See the individual TaskSpaceQuantity's for more information.
*
* <h3>Usage</h3>
*
* We expect you to use this class within a Force::Custom::Implementation, but
* this is not the only option. However, this class can only be used within a
* class that has a realizeTopology method. Here are the necessary steps for
* making use of this class:
*
*  -# Specify your tasks with the TaskSpace::addTask method. If your tasks
*     don't change throughout your simulation, you can do this in a constructor.
*  -# Call TaskSpace::realizeTopology within the owning class' realizeTopology.
*     This is used to allocate cache entries.
*  -# Access the quantities (perhaps in a calcForce method).
*
* Look at the ExampleOperationalSpace.cpp example to see how to use this class.
*
* [1] Khatib, Oussama, et al. "Robotics-based synthesis of human motion."
* Journal of physiology-Paris 103.3 (2009): 211-219.
*/
class SimTK_SIMBODY_EXPORT TaskSpace
{
public:

    //==========================================================================
    // nested classes
    //==========================================================================

    /** An abstract class for common task-space Matrix or Vector quantities.
    *
    * All task-space quantities must be capable of providing their explicit
    * value. After this value is computed, it is cached in the State for
    * efficiency. These classes may optionally provide operators that allow
    * for more efficient computations.
    *
    * The template parameter T is the type of the task-space quantity
    * (e.g., Matrix). The parameter S is used when allocating the cache entry;
    * it is the earliest stage at which the cache entry can be provided.
    */
    template <typename T, Stage::Level S=Stage::Position>
    class SimTK_SIMBODY_EXPORT TaskSpaceQuantity {
    public:
        TaskSpaceQuantity() :
                m_tspace(NULL), m_state(NULL), m_cacheIsValid(false) {
        }

        /** Obtain this quantity explicitly. If possible, use operators
        * instead of this method.
        */
        const T& value() const {
            if (!m_cacheIsValid) {
                updateCache(const_cast<TaskSpaceQuantity*>(this)->m_cache);
                const_cast<TaskSpaceQuantity*>(this)->m_cacheIsValid = true;
            }
            return m_cache;
        }

    protected:

        const State& getState() const {
            return *m_state;
        }

        const TaskSpace* m_tspace;

    private:

        virtual void updateCache(T& cache) const = 0;

        CacheEntryIndex m_cacheIndex;
        bool m_cacheIsValid;
        const State* m_state;
        T m_cache;

        //  Methods that TaskSpace uses.
        friend class TaskSpace;

        void realize(const State& s) {
            m_state = &s;
            m_cacheIsValid = false;
        }

        // The earliest stage at which this quantity can be calculated. This
        // is used for creating lazy cache entries.
        static Stage getEarliestStage() { return S; }
    };

    // Forward declarations.
    class JacobianTranspose;
    class Inertia;
    class DynamicallyConsistentJacobianInverseTranspose;
    class InertiaInverse;
    class NullspaceProjectionTranspose;

    /** Relates task-space velocities to generalized speeds; (nst x nu).
    */
    class SimTK_SIMBODY_EXPORT Jacobian : public TaskSpaceQuantity<Matrix> {
    public:
        Jacobian() {};
        const JacobianTranspose& transpose() const;
        const JacobianTranspose& operator~() { return transpose(); }
        /// Using this operator is likely more efficient than obtaining this
        /// matrix explicitly and performing the multiplication on your own.
        Vector operator*(const Vector& u) const;
        // TODO Matrix_<Vec3> operator*(const Matrix& u) const;
        // TODO Matrix operator*(const Matrix& u) const;
        // TODO Matrix operator*(const NullspaceProjection& N) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** Used to compute task-space forces; (nu x nst).
    */
    class SimTK_SIMBODY_EXPORT JacobianTranspose : public TaskSpaceQuantity<Matrix> {
    public:
        JacobianTranspose() {}
        const Jacobian& transpose() const;
        const Jacobian& operator~() { return transpose(); }
        Vector operator*(const Vector_<Vec3>& f_GP) const;
        Vector operator*(const Vector& f_GP) const;
        Vector operator*(const Vec3& f_GP) const;
        // TODO Matrix operator*(const Matrix_<Vec3>& f_GP) const;
        Matrix operator*(const Matrix& f_GP) const;
        Matrix operator*(const TaskSpace::Inertia& Lambda) const;
        Matrix operator*(
                const TaskSpace::DynamicallyConsistentJacobianInverseTranspose&
                JBarT) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** Task-space inertia matrix; \f$ \Lambda = (J A^{-1} J^T)^{-1} \f$
    * (nst x nst).
    */
    class SimTK_SIMBODY_EXPORT Inertia : public TaskSpaceQuantity<Matrix> {
    public:
        Inertia() {}
        const InertiaInverse& inverse() const;
        Vector operator*(const Vector& a) const;
        Vector operator*(const Vec3& a) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** Inverse of task-space inertia matrix;
    * \f$ \Lambda^{-1} = J M^{-1} J^T \f$ (nst x nst).
    *
    * This is only needed for computing the Inertia matrix.
    */
    class SimTK_SIMBODY_EXPORT InertiaInverse : public TaskSpaceQuantity<Matrix> {
    public:
        InertiaInverse() {}
        const Inertia& inverse() const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** Mass-matrix weighted generalized inverse;
    * \f$ \bar{J} = A^{-1} J^T \Lambda \f$ (nu x nst).
    */
    class SimTK_SIMBODY_EXPORT DynamicallyConsistentJacobianInverse :
        public TaskSpaceQuantity<Matrix> {
    public:
        DynamicallyConsistentJacobianInverse() {}
        const DynamicallyConsistentJacobianInverseTranspose& transpose() const;
        const DynamicallyConsistentJacobianInverseTranspose& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& vec) const;
        Matrix operator*(const Matrix& mat) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** (nst x nu).
    */
    class SimTK_SIMBODY_EXPORT DynamicallyConsistentJacobianInverseTranspose :
        public TaskSpaceQuantity<Matrix> {
    public:
        DynamicallyConsistentJacobianInverseTranspose() {}
        const DynamicallyConsistentJacobianInverse& transpose() const;
        const DynamicallyConsistentJacobianInverse& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& g) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** Includes Coriolis forces and the like;
    * \f$ \mu = \bar{J}^T b - \Lambda \dot{J} u \f$ (nst x 1).
    */
    class SimTK_SIMBODY_EXPORT InertialForces : public TaskSpaceQuantity<Vector, Stage::Velocity> {
    public:
        InertialForces() {}
        Vector operator+(const Vector& f) const;
    private:
        void updateCache(Vector& cache) const OVERRIDE_11;
    };

    /** The task-space forces arising from gravity;
    * \f$ p = \bar{J}^T g \f$ (nst x 1).
    */
    class SimTK_SIMBODY_EXPORT Gravity : public TaskSpaceQuantity<Vector> {
    public:
        Gravity() {}
        Vector operator+(const Vector& f) const;
        /** The joint-space gravity forces (nu x 1).
        */
        Vector systemGravity() const;
        /** A shorthand for the joint-space gravity forces (nu x 1).
        */
        Vector g() const { return systemGravity(); }
    private:
        void updateCache(Vector& cache) const OVERRIDE_11;
    };

    /** Used to prioritize tasks; \f$ N = I - \bar{J} J \f$ (nu x nu).
    */
    class SimTK_SIMBODY_EXPORT NullspaceProjection :
            public TaskSpaceQuantity<Matrix> {
    public:
        NullspaceProjection() {}
        const NullspaceProjectionTranspose& transpose() const;
        const NullspaceProjectionTranspose& operator~() const
        { return transpose(); }
        Vector operator* (const Vector& vec) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    /** (nu x nu).
    */
    class SimTK_SIMBODY_EXPORT NullspaceProjectionTranspose :
            public TaskSpaceQuantity<Matrix> {
    public:
        NullspaceProjectionTranspose() {}
        const NullspaceProjection& transpose() const;
        const NullspaceProjection& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& vec) const;
    private:
        void updateCache(Matrix& cache) const OVERRIDE_11;
    };

    //==========================================================================
    // TaskSpace class
    //==========================================================================

    /**
    * @param[in] matter The matter subsystem being controlled.
    * @param[in] gravity The gravity forces, which should be in the same
    *       System as matter.
    */
    TaskSpace(const SimbodyMatterSubsystem& matter,
              const Force::Gravity&         gravityForce) :
    m_matter(matter), m_gravityForce(gravityForce)
    {}

    /** This \a must be called within realizeTopology of the class that owns
    * this object.
    */
    void realizeTopology(State& state) const {
        const_cast<TaskSpace*>(this)->m_jacobianIndex =
                m_matter.allocateLazyCacheEntry(state,
                        Jacobian::getEarliestStage(),
                        new Value<Jacobian>());
        const_cast<TaskSpace*>(this)->m_jacobianTransposeIndex =
                m_matter.allocateLazyCacheEntry(state,
                        JacobianTranspose::getEarliestStage(),
                        new Value<JacobianTranspose>());
        const_cast<TaskSpace*>(this)->m_inertiaIndex =
                m_matter.allocateLazyCacheEntry(state,
                        Inertia::getEarliestStage(),
                        new Value<Inertia>());
        const_cast<TaskSpace*>(this)->m_inertiaInverseIndex =
                m_matter.allocateLazyCacheEntry(state,
                        InertiaInverse::getEarliestStage(),
                        new Value<InertiaInverse>());
        const_cast<TaskSpace*>(this)->m_jacobianInverseIndex =
                m_matter.allocateLazyCacheEntry(state,
                        DynamicallyConsistentJacobianInverse::getEarliestStage(),
                        new Value<DynamicallyConsistentJacobianInverse>());
        const_cast<TaskSpace*>(this)->m_jacobianInverseTransposeIndex =
                m_matter.allocateLazyCacheEntry(state,
                        DynamicallyConsistentJacobianInverseTranspose::getEarliestStage(),
                        new Value<DynamicallyConsistentJacobianInverseTranspose>());
        const_cast<TaskSpace*>(this)->m_inertialForcesIndex =
                m_matter.allocateLazyCacheEntry(state,
                        InertialForces::getEarliestStage(),
                        new Value<InertialForces>());
        const_cast<TaskSpace*>(this)->m_gravityIndex =
                m_matter.allocateLazyCacheEntry(state,
                        Gravity::getEarliestStage(),
                        new Value<Gravity>());
        const_cast<TaskSpace*>(this)->m_nullspaceIndex =
                m_matter.allocateLazyCacheEntry(state,
                        NullspaceProjection::getEarliestStage(),
                        new Value<NullspaceProjection>());
        const_cast<TaskSpace*>(this)->m_nullspaceTransposeIndex =
                m_matter.allocateLazyCacheEntry(state,
                        NullspaceProjectionTranspose::getEarliestStage(),
                        new Value<NullspaceProjectionTranspose>());
    }

    /** The SimbodyMatterSubsystem that this TaskSpace object uses
    * to compute its quantities.
    */
    const SimbodyMatterSubsystem& getMatterSubsystem() const
    { return m_matter; }

    /// @name Access to TaskSpaceQuantity's.
    /// @{
    TASKSPACEQUANTITY_MEMBERS(Jacobian, jacobian, J);
    TASKSPACEQUANTITY_MEMBERS(JacobianTranspose, jacobianTranspose, JT);
    TASKSPACEQUANTITY_MEMBERS(Inertia, inertia, Lambda);
    TASKSPACEQUANTITY_MEMBERS(InertiaInverse, inertiaInverse, LambdaInv);
    TASKSPACEQUANTITY_MEMBERS(DynamicallyConsistentJacobianInverse, jacobianInverse, JBar);
    TASKSPACEQUANTITY_MEMBERS(DynamicallyConsistentJacobianInverseTranspose, jacobianInverseTranspose, JBarT);
    TASKSPACEQUANTITY_MEMBERS(InertialForces, inertialForces, mu);
    TASKSPACEQUANTITY_MEMBERS(Gravity, gravity, p);
    TASKSPACEQUANTITY_MEMBERS(NullspaceProjection, nullspace, N);
    TASKSPACEQUANTITY_MEMBERS(NullspaceProjectionTranspose, nullspaceTranspose, NT);
    /// @}

    /** Use this to identify Station tasks that have been added to this
    * TaskSpace.
    */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(TaskSpace,StationTaskIndex);

    /// @name Station tasks.
    /// @{
    /** Add a task for a station (point fixed on a rigid body).
    * @param[in] body The index of the body on which the point is fixed.
    * @param[in] station The point of the body, expressed in the body frame.
    */
    void addStationTask(MobilizedBodyIndex body, Vec3 station) {
        m_indices.push_back(body);
        m_stations.push_back(station);
    }

    /** The number of calls to add*Task that have been made (nt).
    */
    unsigned int getNumTasks() const {
        return m_indices.size();
    }

    /** The dimensionality of the task space (nst).
    */
    unsigned int getNumScalarTasks() const {
        return 3 * m_indices.size();
    }
    /// @}

    /// @name Convenience calculations.
    /// @{

    /** Obtain the location and velocity, in the ground frame and expressed
    * in the ground frame, of the station for a given station task.
    */
    void findStationLocationAndVelocityInGround(const State& s,
            StationTaskIndex index,
            const Vec3& stationLocationInBody,
            Vec3& locationInGround, Vec3& velocityInGround) const {
        m_matter.getMobilizedBody(m_indices[index])
                .findStationLocationAndVelocityInGround(s,
                        stationLocationInBody,
                        locationInGround, velocityInGround);
    }

    /** Given accelerations, computes inverse dynamics in the task-space
    * \f$ F = \Lambda F^{*} \mu + p \f$ (nst x 1).
    *
    * This is used to perform feedforward control: given a control law that
    * is specified as task-space accelerations, this provides the task-space
    * forces that achieve those accelerations.
    */
    Vector calcInverseDynamics(const State& s, const Vector& taskAccelerations)
        const;

    /** The joint-space gravity forces (nu x 1).
    */
    Vector g(const State& s) const { return getGravity(s).g(); }
    /// @}

private:

    const Array_<MobilizedBodyIndex>& getMobilizedBodyIndices() const
    { return m_indices; }
    const Array_<Vec3>& getStations() const { return m_stations; }

    //==========================================================================
    // Member variables.
    //==========================================================================

    const SimbodyMatterSubsystem& m_matter;
    const Force::Gravity& m_gravityForce;

    // For station tasks.
    Array_<MobilizedBodyIndex> m_indices;
    Array_<Vec3> m_stations;

};


//==============================================================================
// Jacobian
//==============================================================================
void TaskSpace::Jacobian::updateCache(Matrix& cache) const
{
    m_tspace->getMatterSubsystem().calcStationJacobian(
            getState(),
            m_tspace->getMobilizedBodyIndices(),
            m_tspace->getStations(),
            cache);
}

const TaskSpace::JacobianTranspose& TaskSpace::Jacobian::transpose() const
{
    return m_tspace->getJacobianTranspose(getState());
}

Vector TaskSpace::Jacobian::operator*(const Vector& u) const
{
    Vector_<Vec3> Ju;
    m_tspace->getMatterSubsystem().multiplyByStationJacobian(getState(),
            m_tspace->getMobilizedBodyIndices(), m_tspace->getStations(),
            u, Ju);

    // Convert to a Vector.
    Vector out(3 * Ju.size());
    for (unsigned int i = 0; i < Ju.size(); ++i) {
        out[3 * i] = Ju[i][0];
        out[3 * i + 1] = Ju[i][1];
        out[3 * i + 2] = Ju[i][2];
    }

    return out;
}


//==============================================================================
// JacobianTranspose
//==============================================================================
void TaskSpace::JacobianTranspose::updateCache(Matrix& cache) const
{
    cache = transpose().value().transpose();
}

const TaskSpace::Jacobian& TaskSpace::JacobianTranspose::transpose() const
{
    return m_tspace->getJacobian(getState());
}

Vector TaskSpace::JacobianTranspose::operator*(const Vector_<Vec3>& f_GP) const
{
    Vector f;
    m_tspace->getMatterSubsystem().multiplyByStationJacobianTranspose(
            getState(),
            m_tspace->getMobilizedBodyIndices(),
            m_tspace->getStations(),
            f_GP,
            f);
    return f;
}

Vector TaskSpace::JacobianTranspose::operator*(const Vector& f_GP) const
{
    unsigned int nIn = f_GP.size();
    SimTK_APIARGCHECK1_ALWAYS(nIn % 3 == 0,
            "TaskSpace::JacobianTranspose", "operator*",
            "Length of f_GP, which is %i, is not divisible by 3.", nIn);

    unsigned int nOut = nIn / 3;

    // Create the Vector_<Vec3>.
    // TODO debug, or look for methods that already do this.
    Vector_<Vec3> my_f_GP(nOut);
    for (unsigned int i = 0; i < nOut; ++i)
    {
        // getAs is just a recast; doesn't copy.
        my_f_GP[i] = Vec3::getAs(&f_GP[3 * i]);
    }

    // Perform the multiplication.
    return operator*(my_f_GP);
}

Vector TaskSpace::JacobianTranspose::operator*(const Vec3& f_GP) const
{
   return operator*(Vector_<Vec3>(1, f_GP));
}

Matrix TaskSpace::JacobianTranspose::operator*(const Matrix& f_GP) const
{
    unsigned int nrow = getState().getNU();
    unsigned int ncol = f_GP.ncol();

    Matrix out(nrow, ncol);
    for (unsigned int j = 0; j < ncol; ++j)
    {
        // TODO is this cast inefficient? Is it copying?
        out(j) = operator*(Vector(f_GP(j)));
    }

    return out;
}

Matrix TaskSpace::JacobianTranspose::operator*(
        const TaskSpace::Inertia& Lambda) const
{
    // TOOD could be more efficient.
    return operator*(Lambda.value());
}

Matrix TaskSpace::JacobianTranspose::operator*(
        const TaskSpace::DynamicallyConsistentJacobianInverseTranspose& JBarT) const
{
    return operator*(JBarT.value());
}

//==============================================================================
// Inertia
//==============================================================================
void TaskSpace::Inertia::updateCache(Matrix& cache) const
{
    FactorLU inertiaInverse(m_tspace->getInertiaInverse(getState()).value());
    inertiaInverse.inverse(cache);
}

const TaskSpace::InertiaInverse& TaskSpace::Inertia::inverse() const
{
    return m_tspace->getInertiaInverse(getState());
}

Vector TaskSpace::Inertia::operator*(const Vector& a) const
{
    return value() * a;
}

Vector TaskSpace::Inertia::operator*(const Vec3& a) const
{
    return operator*(Vector(a));
}

//==============================================================================
// InertiaInverse
//==============================================================================
void TaskSpace::InertiaInverse::updateCache(Matrix& cache) const
{
    const SimbodyMatterSubsystem& matter = m_tspace->getMatterSubsystem();

    const JacobianTranspose& JT = m_tspace->JT(getState());
    // TODO const Matrix& JT = m_tspace.JT().value();
    const Matrix& J = m_tspace->J(getState()).value();
    /* TODO 
    // TODO cache the result.

    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix J = m_tspace.getJacobian().value();

    Matrix MInvJt(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        matter.multiplyByMInv(m_tspace.getState(), J.transpose()(j), MInvJt(j));
    }

    updCacheValue() = J * MInvJt;
    */

    unsigned int nt = m_tspace->getNumTasks();
    unsigned int nst = m_tspace->getNumScalarTasks();
    unsigned int nu = getState().getNU();

    Matrix& inertiaInverse = cache;
    inertiaInverse.resize(nst, nst);

    // Create temporary variables.
    Vector Jtcol(nu);
    Vector MInvJtcol(nu);
    Vector_<Vec3> JMInvJt_j(nt);

    // f_GP is used to pluck out one column at a time of Jt. Exactly one
    // element at a time of f_GP will be 1, the rest are 0.
    Vector f_GP(nst, Real(0));

    for (unsigned int j = 0; j < nst; ++j)
    {
        f_GP[j] = 1;
        Jtcol = JT * f_GP;
        f_GP[j] = 0;

        matter.multiplyByMInv(getState(), Jtcol, MInvJtcol);

        // TODO replace with operator.
        inertiaInverse(j) = J * MInvJtcol;
        /* TODO
        matter.multiplyByStationJacobian(m_tspace.getState(),
                m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
                MInvJtcol, JMInvJt_j);

        inertiaInverse(j) = JMInvJt_j;
        */
    }
}

const TaskSpace::Inertia& TaskSpace::InertiaInverse::inverse() const
{
    return m_tspace->getInertia(getState());
}


//==============================================================================
// DynamicallyConsistentJacobianInverse
//==============================================================================
void TaskSpace::DynamicallyConsistentJacobianInverse::updateCache(Matrix& cache)
    const
{
    const JacobianTranspose& JT = m_tspace->getJacobianTranspose(getState());
    const Inertia& Lambda = m_tspace->getInertia(getState());

    // TODO inefficient?
    Matrix JtLambda = JT * Lambda;

    unsigned int nst = m_tspace->getNumScalarTasks();
    unsigned int nu = getState().getNU();

    Matrix& Jbar = cache;
    Jbar.resize(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        m_tspace->getMatterSubsystem().multiplyByMInv(getState(),
                JtLambda(j), Jbar(j));
    }
}

const TaskSpace::DynamicallyConsistentJacobianInverseTranspose&
TaskSpace::DynamicallyConsistentJacobianInverse::transpose() const
{
    return m_tspace->getDynamicallyConsistentJacobianInverseTranspose(getState());
}

Vector TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Vector& vec) const
{
    const JacobianTranspose& JT = m_tspace->getJacobianTranspose(getState());
    const Inertia& Lambda = m_tspace->getInertia(getState());

    // TODO where is this even used? TODO test this.

    Vector JBarvec;
    m_tspace->getMatterSubsystem().multiplyByMInv(getState(),
            JT * (Lambda * vec), JBarvec);
    return  JBarvec;
}

Matrix TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Matrix& mat) const
{
    unsigned int nrow = getState().getNU();
    unsigned int ncol = mat.ncol();

    Matrix out(nrow, ncol);
    for (unsigned int j = 0; j < ncol; ++j)
    {
        out(j) = operator*(mat(j).getAsVector());
    }

    return out;
}


//==============================================================================
// DynamicallyConsistentJacobianInverseTranspose
//==============================================================================
void TaskSpace::DynamicallyConsistentJacobianInverseTranspose::updateCache(
            Matrix& cache) const
{

    cache = transpose().value().transpose();
}

const TaskSpace::DynamicallyConsistentJacobianInverse&
TaskSpace::DynamicallyConsistentJacobianInverseTranspose::transpose() const
{
    return m_tspace->getDynamicallyConsistentJacobianInverse(getState());
}

Vector TaskSpace::DynamicallyConsistentJacobianInverseTranspose::operator*(
        const Vector& g) const
{
    // TODO inefficient. can we have an MInvT operator??
    return value() * g;
}


//==============================================================================
// InertialForces
//==============================================================================
void TaskSpace::InertialForces::updateCache(Vector& cache) const
{
    Vector jointSpaceInertialForces;
    m_tspace->getMatterSubsystem().calcResidualForceIgnoringConstraints(
            getState(), Vector(0), Vector_<SpatialVec>(0), Vector(0),
            jointSpaceInertialForces);

    Vector JDotu;
    m_tspace->getMatterSubsystem().calcBiasForStationJacobian(
            getState(),
            m_tspace->getMobilizedBodyIndices(), m_tspace->getStations(),
            JDotu);

    const DynamicallyConsistentJacobianInverseTranspose& JBarT =
        m_tspace->JBarT(getState());
    const Vector& b = jointSpaceInertialForces;
    const Inertia& Lambda = m_tspace->Lambda(getState());

    cache = JBarT * b - Lambda * JDotu;
}

Vector TaskSpace::InertialForces::operator+(const Vector& f) const
{
    return value() + f;
}

Vector operator+(const Vector& f, const TaskSpace::InertialForces& p)
{
    return f + p.value();
}


//==============================================================================
// Gravity
//==============================================================================
void TaskSpace::Gravity::updateCache(Vector& cache) const
{
    cache = m_tspace->JBarT(getState()) * systemGravity();
}

Vector TaskSpace::Gravity::operator+(const Vector& f) const
{
    return value() + f;
}

// TODO global function.
Vector operator+(const Vector& f, const TaskSpace::Gravity& p)
{
    return f + p.value();
}

Vector TaskSpace::Gravity::systemGravity() const
{
    // TODO where does this go? Make a separate class?
    Vector g;
    m_tspace->getMatterSubsystem().multiplyBySystemJacobianTranspose(
            getState(),
            m_tspace->m_gravityForce.getBodyForces(getState()),
            g);
    // Negate, since we want the 'g' that appears on the same side of the
    // equations of motion as does the mass matrix. That is, M udot + C + g = F
    return -g;
}


//==============================================================================
// NullspaceProjection
//==============================================================================
void TaskSpace::NullspaceProjection::updateCache(Matrix& cache) const
{
    cache = transpose().value().transpose();
}

const TaskSpace::NullspaceProjectionTranspose&
TaskSpace::NullspaceProjection::transpose() const
{
    return m_tspace->getNullspaceProjectionTranspose(getState());
}

Vector TaskSpace::NullspaceProjection::operator*(const Vector& vec)
    const
{
    return vec - (m_tspace->JBar(getState()) * (m_tspace->J(getState()) * vec));
}


//==============================================================================
// NullspaceProjectionTranspose
//==============================================================================
void TaskSpace::NullspaceProjectionTranspose::updateCache(Matrix& cache) const
{
    cache = 1 - (m_tspace->JT(getState()) * m_tspace->JBarT(getState()));
}

const TaskSpace::NullspaceProjection&
TaskSpace::NullspaceProjectionTranspose::transpose() const
{
    return m_tspace->getNullspaceProjection(getState());
}

Vector TaskSpace::NullspaceProjectionTranspose::operator*(const Vector& vec)
    const
{
    return vec - (m_tspace->JT(getState()) * (m_tspace->JBarT(getState()) * vec));
}


//==============================================================================
// TaskSpace
//==============================================================================
// TODO account for applied forces? velocities?
Vector TaskSpace::calcInverseDynamics(const State& s, const Vector& taskAccelerations) const
{
    return Lambda(s) * taskAccelerations + mu(s) + p(s);
}

} // end namespace
