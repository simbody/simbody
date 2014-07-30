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
        // TODO Vector_<Vec3> operator*(const Vector& u) const;
        // TODO Vector operator*(const Vector& u) const;
        // TODO Matrix_<Vec3> operator*(const Matrix& u) const;
        // TODO Matrix operator*(const Matrix& u) const;
    };

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
        // TODO Matrix operator*(const Matrix& f_GP) const;
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

    const Jacobian getJacobian() const { return m_jacobian; }

    const JacobianTranspose getJacobianTranspose() const
    { return m_jacobianTranspose; }

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
    for (unsigned int i = 0; i < nIn; i += 3)
    {
        my_f_GP[i / 3] = Vec3(f_GP[i], f_GP[i + 1], f_GP[i + 2]);
    }

    // Perform the multiplication.
    return operator*(my_f_GP);
}

Vector TaskSpace::JacobianTranspose::operator*(const Vec3& f_GP) const
{
   return operator*(Vector_<Vec3>(1, f_GP));
}

} // end namespace
