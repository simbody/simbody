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
*  -# TODO setState.
*  -# Access the quantities (perhaps in a calcForce method).
*
* TODO examples.
*
*  [1] Khatib, Oussama, et al. "Robotics-based synthesis of human motion."
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
        TaskSpaceQuantity(const TaskSpace& tspace) : m_tspace(tspace) {}
        /** Obtain this quantity explicitly. If possible, use operators
        * instead of this method.
        */
        const T& value() const {
            realizeCacheEntry();
            return getCacheValue();
        }
        void realizeTopology(State& state) const
        {
            const_cast<TaskSpaceQuantity*>(this)->m_cacheIndex =
                m_tspace.getMatterSubsystem().allocateLazyCacheEntry(state,
                        S, new Value<T>());
        }
        virtual void realizeCacheEntry() const = 0;

    protected:
        bool isCacheValueRealized() const {
            return m_tspace.isCacheValueRealized(m_cacheIndex);
        }
    
        void markCacheValueRealized() const {
            return m_tspace.markCacheValueRealized(m_cacheIndex);
        }
    
        const T& getCacheValue() const {
            return Value<T>::downcast(m_tspace.getCacheEntry(m_cacheIndex));
        }
    
        T& updCacheValue() const {
            return Value<T>::updDowncast(m_tspace.updCacheEntry(m_cacheIndex));
        }

        const TaskSpace& m_tspace;

    private:
        CacheEntryIndex m_cacheIndex;
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
        Jacobian(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const JacobianTranspose& transpose() const;
        const JacobianTranspose& operator~() { return transpose(); }
        /// Using this operator is likely efficient than obtaining this matrix
        /// explicitly and performing the multiplication on your own.
        Vector operator*(const Vector& u) const;
        // TODO Matrix_<Vec3> operator*(const Matrix& u) const;
        // TODO Matrix operator*(const Matrix& u) const;
        // TODO Matrix operator*(const NullspaceProjection& N) const;
    };

    /** Used to compute task-space forces; (nu x nst).
    */
    class SimTK_SIMBODY_EXPORT JacobianTranspose : public TaskSpaceQuantity<Matrix> {
    public:
        JacobianTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
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
    };

    /** Task-space inertia matrix; \f$ \Lambda = (J A^{-1} J^T)^{-1} \f$
    * (nst x nst).
    */
    class SimTK_SIMBODY_EXPORT Inertia : public TaskSpaceQuantity<Matrix> {
    public:
        Inertia(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const InertiaInverse& inverse() const;
        // TODO input is a sort of acceleration
        Vector operator*(const Vector& a) const;
        // TODO convenience.
        Vector operator*(const Vec3& a) const;
    };

    /** Inverse of task-space inertia matrix;
    * \f$ \Lambda^{-1} = J M^{-1} J^T \f$ (nst x nst).
    *
    * This is only needed for computing the Inertia matrix.
    */
    class SimTK_SIMBODY_EXPORT InertiaInverse : public TaskSpaceQuantity<Matrix> {
    public:
        InertiaInverse(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const Inertia& inverse() const;
    };

    /** Mass-matrix weighted generalized inverse;
    * \f$ \bar{J} = A^{-1} J^T \Lambda \f$ (nu x nst).
    */
    class SimTK_SIMBODY_EXPORT DynamicallyConsistentJacobianInverse :
        public TaskSpaceQuantity<Matrix> {
    public:
        DynamicallyConsistentJacobianInverse(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const DynamicallyConsistentJacobianInverseTranspose& transpose() const;
        const DynamicallyConsistentJacobianInverseTranspose& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& vec) const;
        Matrix operator*(const Matrix& mat) const;
    };

    /** (nst x nu).
    */
    class SimTK_SIMBODY_EXPORT DynamicallyConsistentJacobianInverseTranspose :
        public TaskSpaceQuantity<Matrix> {
    public:
        DynamicallyConsistentJacobianInverseTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const DynamicallyConsistentJacobianInverse& transpose() const;
        const DynamicallyConsistentJacobianInverse& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& g) const;
    };

    /** Includes Coriolis forces and the like;
    * \f$ \mu = \bar{J}^T b - \Lambda \dot{J} u \f$ (nst x 1).
    */
    class SimTK_SIMBODY_EXPORT InertialForces : public TaskSpaceQuantity<Vector, Stage::Velocity> {
    public:
        InertialForces(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        Vector operator+(const Vector& f) const;
    };

    /** The task-space forces arising from gravity;
    * \f$ p = \bar{J}^T g \f$ (nst x 1).
    */
    class SimTK_SIMBODY_EXPORT Gravity : public TaskSpaceQuantity<Vector> {
    public:
        Gravity(const TaskSpace& tspace, const Force::Gravity& gravity) :
            TaskSpaceQuantity(tspace), m_gravity(gravity) {}
        void realizeCacheEntry() const OVERRIDE_11;
        Vector operator+(const Vector& f) const;
        /** The joint-space gravity forces (nu x 1).
        */
        Vector systemGravity() const;
        /** A shorthand for the joint-space gravity forces (nu x 1).
        */
        Vector g() const { return systemGravity(); }
    private:
        const Force::Gravity& m_gravity;
    };

    /** Used to prioritize tasks; \f$ N = I - \bar{J} J \f$ (nu x nu).
    */
    class SimTK_SIMBODY_EXPORT NullspaceProjection : public TaskSpaceQuantity<Matrix> {
    public:
        NullspaceProjection(const TaskSpace & tspace) :
            TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const NullspaceProjectionTranspose& transpose() const;
        const NullspaceProjectionTranspose& operator~() const
        { return transpose(); }
        Vector operator* (const Vector& vec) const;
    };

    /** (nu x nu).
    */
    class SimTK_SIMBODY_EXPORT NullspaceProjectionTranspose : public TaskSpaceQuantity<Matrix> {
    public:
        NullspaceProjectionTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        void realizeCacheEntry() const OVERRIDE_11;
        const NullspaceProjection& transpose() const;
        const NullspaceProjection& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& vec) const;
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
              const Force::Gravity&         gravity) :
    m_matter(matter), m_jacobian(*this), m_jacobianTranspose(*this),
    m_inertia(*this), m_inertiaInverse(*this),
    m_jacobianInverse(*this), m_jacobianInverseTranspose(*this),
    m_inertialForces(*this), m_gravity(*this, gravity),
    m_nullspace(*this), m_nullspaceTranspose(*this)
    {}

    /** This \a must be called within realizeTopology of the class that owns
    * this object.
    */
    void realizeTopology(State& state) const {
        m_jacobian.realizeTopology(state);
        m_jacobianTranspose.realizeTopology(state);
        m_inertia.realizeTopology(state);
        m_inertiaInverse.realizeTopology(state);
        m_jacobianInverse.realizeTopology(state);
        m_jacobianInverseTranspose.realizeTopology(state);
        m_gravity.realizeTopology(state);
        m_nullspace.realizeTopology(state);
        m_nullspaceTranspose.realizeTopology(state);
    }

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

    /** TODO
    */
    void setState(const State& state) const {
        // TODO is this okay?
        const_cast<TaskSpace*>(this)->m_state = &state;
    }

    /** TODO
    */
    const State& getState() const {
        if (!m_state)
        {
            // TODO
            throw Exception::Base("State is null.");
        }
        return *m_state;
    }

    const SimbodyMatterSubsystem& getMatterSubsystem() const
    { return m_matter; }
    const Array_<MobilizedBodyIndex>& getMobilizedBodyIndices() const
    { return m_indices; }
    const Array_<Vec3>& getStations() const { return m_stations; }

    /// @name Shorthand access to the task space quantities.
    /// @{
    const Jacobian& J() const { return getJacobian(); }
    const JacobianTranspose& JT() const
    { return getJacobianTranspose(); }
    const Inertia& Lambda() const { return getInertia(); }
    const InertiaInverse& LambdaInv() const { return getInertiaInverse(); }
    const DynamicallyConsistentJacobianInverse& JBar() const
    { return getDynamicallyConsistentJacobianInverse(); }
    const DynamicallyConsistentJacobianInverseTranspose& JBarT() const
    { return getDynamicallyConsistentJacobianInverseTranspose(); }
    const InertialForces mu() const { return getInertialForces(); }
    const Gravity p() const { return getGravity(); }
    const NullspaceProjection& N() const { return getNullspaceProjection(); }
    const NullspaceProjectionTranspose NT() const
    { return getNullspaceProjectionTranspose(); }
    /// @}

    /// @name Access to the task space quantities.
    /// @{
    const Jacobian& getJacobian() const { return m_jacobian; }
    const JacobianTranspose& getJacobianTranspose() const
    { return m_jacobianTranspose; }
    const Inertia& getInertia() const { return m_inertia; }
    const InertiaInverse& getInertiaInverse() const { return m_inertiaInverse; }
    const DynamicallyConsistentJacobianInverse&
        getDynamicallyConsistentJacobianInverse() const
    { return m_jacobianInverse; }
    const DynamicallyConsistentJacobianInverseTranspose&
        getDynamicallyConsistentJacobianInverseTranspose() const
    { return m_jacobianInverseTranspose; }
    const InertialForces& getInertialForces() const { return m_inertialForces; }
    const Gravity& getGravity() const { return m_gravity; }
    const NullspaceProjection& getNullspaceProjection() const
    { return m_nullspace; }
    const NullspaceProjectionTranspose& getNullspaceProjectionTranspose() const
    { return m_nullspaceTranspose; }
    /// @}


    /** Given accelerations, computes inverse dynamics in the task-space
    * \f$ F = \Lambda F^{*} \mu + p \f$ (nst x 1).
    *
    * This is used to perform feedforward control: given a control law that
    * is specified as task-space accelerations, this provides the task-space
    * forces that achieve those accelerations.
    *
    * TODO does not use mu yet.
    */
    Vector calcInverseDynamics(const Vector& taskAccelerations) const;

    /** The joint-space gravity forces (nu x 1).
    */
    Vector g() const { return getGravity().g(); }

private:

    //==========================================================================
    // Cache value helpers.
    //==========================================================================

    bool isCacheValueRealized(CacheEntryIndex index) const {
        return m_matter.isCacheValueRealized(*m_state, index);
    }

    void markCacheValueRealized(CacheEntryIndex index) const {
        return m_matter.markCacheValueRealized(*m_state, index);
    }

    const AbstractValue& getCacheEntry(CacheEntryIndex index) const {
        return m_matter.getCacheEntry(*m_state, index);
    }

    AbstractValue& updCacheEntry(CacheEntryIndex index) const {
        return m_matter.updCacheEntry(*m_state, index);
    }

    //==========================================================================
    // Member variables.
    //==========================================================================

    const SimbodyMatterSubsystem& m_matter;
    const State* m_state;

    Array_<MobilizedBodyIndex> m_indices;
    Array_<Vec3> m_stations;

    Jacobian m_jacobian;
    JacobianTranspose m_jacobianTranspose;
    Inertia m_inertia;
    InertiaInverse m_inertiaInverse;
    DynamicallyConsistentJacobianInverse m_jacobianInverse;
    DynamicallyConsistentJacobianInverseTranspose m_jacobianInverseTranspose;
    InertialForces m_inertialForces;
    Gravity m_gravity;
    NullspaceProjection m_nullspace;
    NullspaceProjectionTranspose m_nullspaceTranspose;

};


//==============================================================================
// Jacobian
//==============================================================================
void TaskSpace::Jacobian::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    m_tspace.getMatterSubsystem().calcStationJacobian(
            m_tspace.getState(),
            m_tspace.getMobilizedBodyIndices(),
            m_tspace.getStations(),
            updCacheValue());

    markCacheValueRealized();
}

const TaskSpace::JacobianTranspose& TaskSpace::Jacobian::transpose() const
{
    return m_tspace.getJacobianTranspose();
}

Vector TaskSpace::Jacobian::operator*(const Vector& u) const
{
    Vector_<Vec3> Ju;
    m_tspace.getMatterSubsystem().multiplyByStationJacobian(m_tspace.getState(),
            m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
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
void TaskSpace::JacobianTranspose::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    // TODO is there any cost to transposing each time?
    updCacheValue() = transpose().value().transpose();

    markCacheValueRealized();
}

const TaskSpace::Jacobian& TaskSpace::JacobianTranspose::transpose() const
{
    return m_tspace.getJacobian();
}

Vector TaskSpace::JacobianTranspose::operator*(const Vector_<Vec3>& f_GP) const
{
    Vector f;
    m_tspace.getMatterSubsystem().multiplyByStationJacobianTranspose(
            m_tspace.getState(),
            m_tspace.getMobilizedBodyIndices(),
            m_tspace.getStations(),
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
    unsigned int nrow = m_tspace.getState().getNU();
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
void TaskSpace::Inertia::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    FactorLU inertiaInverse(m_tspace.getInertiaInverse().value());
    inertiaInverse.inverse(updCacheValue());

    markCacheValueRealized();
}

const TaskSpace::InertiaInverse& TaskSpace::Inertia::inverse() const
{
    return m_tspace.getInertiaInverse();
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
void TaskSpace::InertiaInverse::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    const SimbodyMatterSubsystem& matter = m_tspace.getMatterSubsystem();

    const JacobianTranspose& JT = m_tspace.JT();
    // TODO const Matrix& JT = m_tspace.JT().value();
    const Matrix& J = m_tspace.J().value();
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

    unsigned int nt = m_tspace.getNumTasks();
    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix& inertiaInverse = updCacheValue();
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

        matter.multiplyByMInv(m_tspace.getState(), Jtcol, MInvJtcol);

        // TODO replace with operator.
        inertiaInverse(j) = J * MInvJtcol;
        /* TODO
        matter.multiplyByStationJacobian(m_tspace.getState(),
                m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
                MInvJtcol, JMInvJt_j);

        inertiaInverse(j) = JMInvJt_j;
        */
    }

    markCacheValueRealized();
}

const TaskSpace::Inertia& TaskSpace::InertiaInverse::inverse() const
{
    return m_tspace.getInertia();
}


//==============================================================================
// DynamicallyConsistentJacobianInverse
//==============================================================================
void TaskSpace::DynamicallyConsistentJacobianInverse::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    const JacobianTranspose& JT = m_tspace.getJacobianTranspose();
    const Inertia& Lambda = m_tspace.getInertia();

    // TODO inefficient?
    Matrix JtLambda = JT * Lambda;

    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix& Jbar = updCacheValue();
    Jbar.resize(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        m_tspace.getMatterSubsystem().multiplyByMInv(m_tspace.getState(),
                JtLambda(j), Jbar(j));
    }

    markCacheValueRealized();
}

const TaskSpace::DynamicallyConsistentJacobianInverseTranspose&
TaskSpace::DynamicallyConsistentJacobianInverse::transpose() const
{
    return m_tspace.getDynamicallyConsistentJacobianInverseTranspose();
}

Vector TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Vector& vec) const
{
    const JacobianTranspose& JT = m_tspace.getJacobianTranspose();
    const Inertia& Lambda = m_tspace.getInertia();

    // TODO where is this even used? TODO test this.

    Vector JBarvec;
    m_tspace.getMatterSubsystem().multiplyByMInv(m_tspace.getState(),
            JT * (Lambda * vec), JBarvec);
    return  JBarvec;
}

Matrix TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Matrix& mat) const
{
    unsigned int nrow = m_tspace.getState().getNU();
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
void
TaskSpace::DynamicallyConsistentJacobianInverseTranspose::realizeCacheEntry()
    const
{
    if (isCacheValueRealized())
        return;

    updCacheValue() = transpose().value().transpose();

    markCacheValueRealized();
}

const TaskSpace::DynamicallyConsistentJacobianInverse&
TaskSpace::DynamicallyConsistentJacobianInverseTranspose::transpose() const
{
    return m_tspace.getDynamicallyConsistentJacobianInverse();
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
void TaskSpace::InertialForces::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    Vector jointSpaceInertialForces;
    m_tspace.getMatterSubsystem().calcResidualForceIgnoringConstraints(
            m_tspace.getState(), Vector(0), Vector_<SpatialVec>(0), Vector(0),
            jointSpaceInertialForces);

    Vector JDotu;
    m_tspace.getMatterSubsystem().calcBiasForStationJacobian(
            m_tspace.getState(), 
            m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
            JDotu);

    const DynamicallyConsistentJacobianInverseTranspose& JBarT =
        m_tspace.JBarT();
    const Vector& b = jointSpaceInertialForces;
    const Inertia& Lambda = m_tspace.Lambda();

    updCacheValue() = JBarT * b - Lambda * JDotu;

    markCacheValueRealized();
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
void TaskSpace::Gravity::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    updCacheValue() = m_tspace.JBarT() * systemGravity();

    markCacheValueRealized();
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
    m_tspace.getMatterSubsystem().multiplyBySystemJacobianTranspose(
            m_tspace.getState(),
            m_gravity.getBodyForces(m_tspace.getState()),
            g);
    // Negate, since we want the 'g' that appears on the same side of the
    // equations of motion as does the mass matrix. That is, M udot + C + g = F
    return -g;
}


//==============================================================================
// NullspaceProjection
//==============================================================================
void TaskSpace::NullspaceProjection::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    updCacheValue() = transpose().value().transpose();

    markCacheValueRealized();
}

const TaskSpace::NullspaceProjectionTranspose&
TaskSpace::NullspaceProjection::transpose() const
{
    return m_tspace.getNullspaceProjectionTranspose();
}

Vector TaskSpace::NullspaceProjection::operator*(const Vector& vec)
    const
{
    return vec - (m_tspace.JBar() * (m_tspace.J() * vec));
}


//==============================================================================
// NullspaceProjectionTranspose
//==============================================================================
void TaskSpace::NullspaceProjectionTranspose::realizeCacheEntry() const
{
    if (isCacheValueRealized())
        return;

    updCacheValue() = 1 - (m_tspace.JT() * m_tspace.JBarT());

    markCacheValueRealized();
}

const TaskSpace::NullspaceProjection&
TaskSpace::NullspaceProjectionTranspose::transpose() const
{
    return m_tspace.getNullspaceProjection();
}

Vector TaskSpace::NullspaceProjectionTranspose::operator*(const Vector& vec)
    const
{
    return vec - (m_tspace.JT() * (m_tspace.JBarT() * vec));
}


//==============================================================================
// TaskSpace
//==============================================================================
// TODO account for applied forces? velocities?
Vector TaskSpace::calcInverseDynamics(const Vector& taskAccelerations) const
{
    return Lambda() * taskAccelerations + /* TODO mu() + */ p();
}

} // end namespace
