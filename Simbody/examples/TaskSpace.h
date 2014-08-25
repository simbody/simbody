#include <Simbody.h>

namespace SimTK {

class TaskSpace
{
public:

    //==========================================================================
    // nested classes
    //==========================================================================

    template <typename T>
    class TaskSpaceQuantity
    {
    public:
        TaskSpaceQuantity(const TaskSpace& tspace) : m_tspace(tspace) {}
        virtual T value() const = 0;
    protected:
        const TaskSpace& m_tspace;
    };

    // Forward declaration.
    class JacobianTranspose;
    class Jacobian : public TaskSpaceQuantity<Matrix> {
    public:
        Jacobian(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const JacobianTranspose& transpose() const;
        const JacobianTranspose& operator~() { return transpose(); }
        // TODO Vector_<Vec3> operator*(const Vector& u) const;
        // TODO Vector operator*(const Vector& u) const;
        // TODO Matrix_<Vec3> operator*(const Matrix& u) const;
        // TODO Matrix operator*(const Matrix& u) const;
        // TODO Matrix operator*(const NullspaceProjection& N) const;
    };

    // Forward declaration.
    class Inertia;
    class DynamicallyConsistentJacobianInverseTranspose;
    class JacobianTranspose : public TaskSpaceQuantity<Matrix> {
    public:
        JacobianTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const Jacobian& transpose() const;
        const Jacobian& operator~() { return transpose(); }
        // TODO I want the input to be a Vector.
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
    
    // Forward declaration.
    class InertiaInverse;
    class Inertia : public TaskSpaceQuantity<Matrix> {
    public:
        Inertia(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const InertiaInverse& inverse() const;
        // TODO input is a sort of acceleration
        Vector operator*(const Vector& a) const;
    };

    /// J M^{-1} J^T
    class InertiaInverse : public TaskSpaceQuantity<Matrix> {
    public:
        InertiaInverse(const TaskSpace& tspace) : TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const Inertia& inverse() const;
    };

    /// M^{-1} J^T Inertia
    class DynamicallyConsistentJacobianInverse :
        public TaskSpaceQuantity<Matrix> {
    public:
        DynamicallyConsistentJacobianInverse(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const DynamicallyConsistentJacobianInverseTranspose& transpose() const;
        const DynamicallyConsistentJacobianInverseTranspose& operator~() const
        { return transpose(); }
        Matrix operator*(const Matrix& mat) const;
    };

    class DynamicallyConsistentJacobianInverseTranspose :
        public TaskSpaceQuantity<Matrix> {
    public:
        DynamicallyConsistentJacobianInverseTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const DynamicallyConsistentJacobianInverse& transpose() const;
        const DynamicallyConsistentJacobianInverse& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& g) const;
    };

    class Gravity : public TaskSpaceQuantity<Vector> {
    public:
        Gravity(const TaskSpace& tspace, const Force::Gravity& gravity) :
            TaskSpaceQuantity(tspace), m_gravity(gravity) {}
        Vector value() const OVERRIDE_11;
        Vector operator+(const Vector& f) const;
        Vector systemGravity() const;
        Vector g() const { return systemGravity(); }
    private:
        const Force::Gravity& m_gravity;
    };

    // Forward declaration.
    class NullspaceProjectionTranspose;
    class NullspaceProjection : public TaskSpaceQuantity<Matrix> {
    public:
        NullspaceProjection(const TaskSpace & tspace) :
            TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const NullspaceProjectionTranspose& transpose() const;
        const NullspaceProjectionTranspose& operator~() const
        { return transpose(); }
    };

    class NullspaceProjectionTranspose : public TaskSpaceQuantity<Matrix> {
    public:
        NullspaceProjectionTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const NullspaceProjection& transpose() const;
        const NullspaceProjection& operator~() const
        { return transpose(); }
        Vector operator*(const Vector& vec) const;
    };

    /// TODO class CoriolisForce; requires calcuating b. calcResidual minus
    // gravity.

    //==========================================================================
    // TaskSpace class
    //==========================================================================

    TaskSpace(const SimbodyMatterSubsystem& matter,
              const Force::Gravity&         gravity) :
    m_matter(matter), m_jacobian(*this), m_jacobianTranspose(*this),
    m_inertia(*this), m_inertiaInverse(*this),
    m_jacobianInverse(*this), m_jacobianInverseTranspose(*this),
    m_gravity(*this, gravity), m_nullspace(*this), m_nullspaceTranspose(*this)
    {}

    void addTask(MobilizedBodyIndex body, Vec3 station)
    {
        m_indices.push_back(body);
        m_stations.push_back(station);
    }

    unsigned int getNumTasks() const
    {
        return m_indices.size();
    }

    unsigned int getNumScalarTasks() const
    {
        return 3 * m_indices.size();
    }

    void setState(const State& state) const
    {
        const_cast<TaskSpace*>(this)->m_state = &state;
    }

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
    const Gravity& getGravity() const { return m_gravity; }
    const NullspaceProjection& getNullspaceProjection() const
    { return m_nullspace; }
    const NullspaceProjectionTranspose& getNullspaceProjectionTranspose() const
    { return m_nullspaceTranspose; }

    // shorthands

    const Jacobian& J() const { return getJacobian(); }
    const JacobianTranspose& JT() const
    { return getJacobianTranspose(); }
    const Inertia& Lambda() const { return getInertia(); }
    const InertiaInverse& LambdaInv() const { return getInertiaInverse(); }
    const DynamicallyConsistentJacobianInverse& JBar() const
    { return getDynamicallyConsistentJacobianInverse(); }
    const DynamicallyConsistentJacobianInverseTranspose& JBarT() const
    { return getDynamicallyConsistentJacobianInverseTranspose(); }
    const Gravity p() const { return getGravity(); }
    const NullspaceProjection& N() const { return getNullspaceProjection(); }
    const NullspaceProjectionTranspose NT() const
    { return getNullspaceProjectionTranspose(); }

    // task-space F = Lambda F* + mu + p
    // TODO does not use mu yet.
    Vector calcInverseDynamics(const Vector& taskAccelerations) const;
    Vector g() const { return getGravity().g(); }

private:

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
    Gravity m_gravity;
    NullspaceProjection m_nullspace;
    NullspaceProjectionTranspose m_nullspaceTranspose;

};


//==============================================================================
// Jacobian
//==============================================================================
Matrix TaskSpace::Jacobian::value() const
{
    Matrix JS;
    m_tspace.getMatterSubsystem().calcStationJacobian(
            m_tspace.getState(),
            m_tspace.getMobilizedBodyIndices(),
            m_tspace.getStations(),
            JS);
    return JS;
}

const TaskSpace::JacobianTranspose& TaskSpace::Jacobian::transpose() const
{
    return m_tspace.getJacobianTranspose();
}


//==============================================================================
// JacobianTranspose
//==============================================================================
Matrix TaskSpace::JacobianTranspose::value() const
{
    return transpose().value().transpose();
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
Matrix TaskSpace::Inertia::value() const
{
    FactorLU inertiaInverse(m_tspace.getInertiaInverse().value());
    Matrix inertia;
    inertiaInverse.inverse(inertia);
    return inertia;
}

const TaskSpace::InertiaInverse& TaskSpace::Inertia::inverse() const
{
    return m_tspace.getInertiaInverse();
}

Vector TaskSpace::Inertia::operator*(const Vector& a) const
{
    return value() * a;
}

//==============================================================================
// InertiaInverse
//==============================================================================
Matrix TaskSpace::InertiaInverse::value() const
{
    // TODO cache the result.
    const SimbodyMatterSubsystem& matter = m_tspace.getMatterSubsystem();

    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix J = m_tspace.getJacobian().value();

    Matrix MInvJt(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        matter.multiplyByMInv(m_tspace.getState(), J.transpose()(j), MInvJt(j));
    }

    return J * MInvJt;

    /*
    unsigned int nt = m_tspace.getNumTasks();
    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix inertiaInverse(nst, nst);

    // Create temporary variables.
    Vector Jtcol(nu);
    Vector MInvJtcol(nu);
    Vector_<Vec3> JMInvJt_j(nt);

    // f_GP is used to pluck out one column at a time of Jt. Exactly one
    // element at a time of f_GP will be 1, the rest are 0.
    Vector_<Vec3> f_GP(nt, Vec3(0));

    for (unsigned int j = 0; j < nst; ++j)
    {
        unsigned int iTask = floor(j / 3);
        unsigned iDim = j % 3;

        f_GP[iTask][iDim] = 1;
        matter.multiplyByStationJacobianTranspose( m_tspace.getState(),
                m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
                f_GP, Jtcol);
        f_GP[iTask][iDim] = 0;

        matter.multiplyByMInv(m_tspace.getState(), Jtcol, MInvJtcol);

        matter.multiplyByStationJacobian(m_tspace.getState(),
                m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
                MInvJtcol, JMInvJt_j);
        inertiaInverse(j) = JMInvJt_j;
    }

    return inertiaInverse;
    */
}

const TaskSpace::Inertia& TaskSpace::InertiaInverse::inverse() const
{
    return m_tspace.getInertia();
}


//==============================================================================
// DynamicallyConsistentJacobianInverse
//==============================================================================
Matrix TaskSpace::DynamicallyConsistentJacobianInverse::value() const
{
    const JacobianTranspose& JT = m_tspace.getJacobianTranspose();
    const Inertia& Lambda = m_tspace.getInertia();

    // TODO inefficient?
    Matrix JtLambda = JT * Lambda;

    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix Jbar(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        m_tspace.getMatterSubsystem().multiplyByMInv(m_tspace.getState(),
                JtLambda(j), Jbar(j));
    }

    return Jbar;
}

const TaskSpace::DynamicallyConsistentJacobianInverseTranspose&
TaskSpace::DynamicallyConsistentJacobianInverse::transpose() const
{
    return m_tspace.getDynamicallyConsistentJacobianInverseTranspose();
}

Matrix TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Matrix& mat) const
{
    return value() * mat;
}


//==============================================================================
// DynamicallyConsistentJacobianInverseTranspose
//==============================================================================
Matrix TaskSpace::DynamicallyConsistentJacobianInverseTranspose::value() const
{
    return transpose().value().transpose();
}

const TaskSpace::DynamicallyConsistentJacobianInverse&
TaskSpace::DynamicallyConsistentJacobianInverseTranspose::transpose() const
{
    return m_tspace.getDynamicallyConsistentJacobianInverse();
}

Vector TaskSpace::DynamicallyConsistentJacobianInverseTranspose::operator*(
        const Vector& g) const
{
    // TODO inefficient
    return value() * g;
}


//==============================================================================
// Gravity
//==============================================================================
Vector TaskSpace::Gravity::value() const
{

    return m_tspace.JBarT() * systemGravity();
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
// Gravity
//==============================================================================
// TODO account for applied forces? velocities?
Vector TaskSpace::calcInverseDynamics(const Vector& taskAccelerations) const
{
    return Lambda() * taskAccelerations + p();
}


//==============================================================================
// NullspaceProjection
//==============================================================================
Matrix TaskSpace::NullspaceProjection::value() const
{
    return transpose().value().transpose();
}

const TaskSpace::NullspaceProjectionTranspose&
TaskSpace::NullspaceProjection::transpose() const
{
    return m_tspace.getNullspaceProjectionTranspose();
}


//==============================================================================
// NullspaceProjectionTranspose
//==============================================================================
Matrix TaskSpace::NullspaceProjectionTranspose::value() const
{
    return -(m_tspace.JT() * m_tspace.JBarT()) + 1;
}

const TaskSpace::NullspaceProjection&
TaskSpace::NullspaceProjectionTranspose::transpose() const
{
    return m_tspace.getNullspaceProjection();
}

Vector TaskSpace::NullspaceProjectionTranspose::operator*(const Vector& vec)
    const
{
    // TODO
    return value() * vec;
}

} // end namespace
