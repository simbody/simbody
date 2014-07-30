#include <Simbody.h>

namespace SimTK {

class TaskSpace
{
public:

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
    };

    class JacobianTranspose : public TaskSpaceQuantity<Matrix> {
    public:
        JacobianTranspose(const TaskSpace& tspace) :
            TaskSpaceQuantity(tspace) {}
        Matrix value() const OVERRIDE_11;
        const Jacobian& transpose() const;
        const Jacobian& operator~() { return transpose(); }
        // TODO I want the input to be a Vector.
        Vector operator*(const Vector_<Vec3>& f_GP);
    };
    // TODO class Inertia;
    // TODO // TODO class CoriolisForce;
    // TODO class GravityForce;
    // TODO class NullspaceProjectionTranspose;

    TaskSpace(const SimbodyMatterSubsystem& matter) :
    m_matter(matter), m_jacobian(*this), m_jacobianTranspose(*this)
    // TODO m_inertia(*this), m_gravity(*this), m_nullspace(*this)
    {}

    void addTask(MobilizedBodyIndex body, Vec3 station)
    {
        m_indices.push_back(body);
        m_stations.push_back(station);
    }

    void setState(const State& state) { m_state = &state; }

    const State& getState() const { return *m_state; }

    const SimbodyMatterSubsystem& getMatterSubsystem() const
    { return m_matter; }

    const Array_<MobilizedBodyIndex>& getMobilizedBodyIndices() const
    { return m_indices; }

    const Array_<Vec3>& getStations() const { return m_stations; }

private:

    const SimbodyMatterSubsystem& m_matter;
    const State* m_state;

    Array_<MobilizedBodyIndex> m_indices;
    Array_<Vec3> m_stations;

    Jacobian m_jacobian;
    JacobianTranspose m_jacobianTranspose;
    // TODO Inertia m_inertia;
    // TODO GravityForce m_gravity;
    // TODO NullspaceProjectionTranspose m_nullspace;

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
    return m_tspace.m_jacobianTranspose;
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
    return m_tspace.m_jacobian;
}

Vector TaskSpace::JacobianTranspose::operator*(const Vector_<Vec3>& f_GP)
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

} // end namespace
