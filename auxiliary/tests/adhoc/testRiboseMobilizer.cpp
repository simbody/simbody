#include "SimTKsimbody_aux.h"
#include <iostream>

using namespace std;
using namespace SimTK;

// eliding boost::units expressions to make a more portable example --cmb
// pay no attention to these types and units.  They are hollow echo of
// the use of boost::units in an earlier version of this example
typedef Real angle_t;
typedef Real length_t;
typedef Real mass_t;
typedef Vec3 location_t;
Real degrees = 3.14159/180.0;
Real radians = 1.0;
Real angstroms = 0.10;
Real nanometers = 1.0;
Real daltons = 1.0;

class ZeroFunction : public Function<1>::Constant {
public:
        ZeroFunction() : Function<1>::Constant(Vec1(0), 0) {}
};

class UnityFunction : public Function<1>::Constant {
public:
    UnityFunction() : Function<1>::Constant(Vec1(1.0), 0) {}
};

class IdentityFunction : public Function<1> {
public:
    IdentityFunction() {}

	Vec<1> calcValue(const Vector& x) const{
        return Vec<1>(x[0]);
    }

	Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
        if (derivComponents.size() == 1) return Vec<1>(1.0);
        else return Vec<1>(0.0);
    }

	int getArgumentSize() const {return 1;}
	
	int getMaxDerivativeOrder() const {return 1000;}
};

// Difference between two components of a two component vector
// Used as a target function for equality coupling constraint
class DifferenceFunction : public Function<1>::Linear {
public:
    DifferenceFunction() : Function<1>::Linear( createCoefficients() ) {}

private:
    static Vector_< Vec<1> > createCoefficients() {
        Vector_< Vec<1> > answer(3);
        answer[0] = Vec<1>(1.0); // ax
        answer[1] = Vec<1>(-1.0); // -by
        answer[2] = Vec<1>(0.0); // +c
        return answer;
    }
};


/// Function that when zero, ensures that three variables are equal
/// f(x,y,z) = (x-y)^2 + (x-z)^2 + (y-z)^2
/// df/dx = 2(x-y) + 2(x-z)
/// df/dx^2 = 4
/// df/dxdy = -2
class ThreeDifferencesFunction : public Function<1> {
public:
    Vec<1> calcValue(const Vector& x) const 
    {
        assert( 3 == x.size() );

        Real dxy(x[0] - x[1]);
        Real dxz(x[0] - x[2]);
        Real dyz(x[1] - x[2]);

        return Vec<1>(dxy*dxy + dxz*dxz + dyz*dyz);
    }

	Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const
    {
		Vec<1> deriv(0);
	
		assert(3 == x.size());
		
        int derivOrder = derivComponents.size();
        assert(1 <= derivOrder);

        if (derivOrder == 1) 
        {
            // too clever
            //    df/dx
            //  = 2(x-y) + 2(x-z) 
            //  = 2x - 2y + 2x - 2z
            //  = 4x - 2y - 2z
            //  = 6x - 2x - 2y - 2z
            deriv = 2 * (3 * x[derivComponents[0]] -x[0] -x[1] -x[2]);
        }
        else if (derivOrder == 2) 
        {
            if (derivComponents[0] == derivComponents[1])
                deriv = 4.0; // df/dx^2
            else
                deriv = -2.0; // df/dxdy
        } 
        else ; // all derivatives higher than two are zero

        return deriv;
    }

	int getArgumentSize() const{
		return 3;
	}
	
	int getMaxDerivativeOrder() const{
		return 1000;
	}
};

/// Implements a simple functional relationship, y = amplitude * sin(x - phase)
class SinusoidFunction : public Function<1> {
private:
	angle_t amplitude;
	angle_t phase;
public:

	//Default constructor
	SinusoidFunction()
		: amplitude(180.0*degrees), phase(0.0*degrees) {}
	
	//Convenience constructor to specify the slope and Y-intercept of the linear r
	SinusoidFunction(angle_t amp, angle_t phi)
	: amplitude(amp), phase(phi) {}
	
	Vec<1> calcValue(const Vector& x) const{
        assert( 1 == x.size() );
        return Vec<1>( angle_t(amplitude*sin(x[0]*radians - phase)) );
	}
	
	Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
		Vec<1> deriv(0);
	
		assert(1 == x.size());
		
        // Derivatives repeat after 4
        int derivOrder = derivComponents.size() % 4;

		// Derivatives 1, 5, 9, 13, ... are cos()
		if      ( 1 == derivOrder ) {
			deriv[0] = angle_t(amplitude*cos(x[0]*radians - phase));
		}
		// Derivatives 2, 6, 10, 14, ... are -sin()
		else if ( 2 == derivOrder ) {
			deriv[0] = angle_t(-amplitude*sin(x[0]*radians - phase));
		}
		// Derivatives 3, 7, 11, 15, ... are -cos()
		else if ( 3 == derivOrder ) {
			deriv[0] = angle_t(-amplitude*cos(x[0]*radians - phase));
		}
		// Derivatives 0, 4, 8, 12, ... are sin()
		else if ( 0 == derivOrder ) {
			deriv[0] = angle_t(amplitude*sin(x[0]*radians - phase));
		}
		else assert(false);
		
		return deriv;
	}
	
	int getArgumentSize() const{
		return 1;
	}
	
	int getMaxDerivativeOrder() const{
		return 1000;
	}
};

// The pupose of TestPinMobilizer is to prove that I have understood the
// basic syntax of Function-based mobilizers
class TestPinMobilizer : public MobilizedBody::FunctionBased
{
public:
        TestPinMobilizer(
                        MobilizedBody& parent,
                        const Transform& inbFrame,
                        const Body& body,
                        const Transform& outbFrame
                        )
        : MobilizedBody::FunctionBased(
                        parent,
                        inbFrame,
                        body,
                        outbFrame,
                        1,
                        createFunctions(),
                        createCoordIndices()
                        )
        {}

private:
    static std::vector< const Function<1>* > createFunctions() {
        std::vector< const Function<1>* > functions;
        functions.push_back(new ZeroFunction ); // x rotation
        functions.push_back(new ZeroFunction ); // y rotation

        // Z-rotation is the only degree of freedom for this mobilizer
        functions.push_back(new IdentityFunction ); // z rotation

        functions.push_back(new ZeroFunction ); // x translation
        functions.push_back(new ZeroFunction ); // y translation
        functions.push_back(new ZeroFunction ); // z translation

        return functions;        
    }

    static std::vector< std::vector<int> > createCoordIndices() {
        std::vector<int> oneCoordinateVec;
        oneCoordinateVec.push_back(0); // first and only generalized coordinate index

        std::vector< std::vector<int> > coordIndices;
        // All six functions take the one generalized coordinate
        coordIndices.push_back(std::vector<int>()); // 1, empty
        coordIndices.push_back(std::vector<int>()); // 2, empty
        coordIndices.push_back(oneCoordinateVec);   // 3, one coordinate, zero
        coordIndices.push_back(std::vector<int>()); // 4, empty
        coordIndices.push_back(std::vector<int>()); // 5, empty
        coordIndices.push_back(std::vector<int>()); // 6, empty

        return coordIndices;
    }
};

class PseudorotationMobilizer : public MobilizedBody::FunctionBased
{
public:
        PseudorotationMobilizer(
                        MobilizedBody& parent,
                        const Transform& inbFrame,
                        const Body& body,
                        const Transform& outbFrame,
                        angle_t amplitude,
                        angle_t phase
                        )
        : MobilizedBody::FunctionBased(
                        parent,
                        inbFrame,
                        body,
                        outbFrame,
                        1,
                        createFunctions(amplitude, phase),
                        createCoordIndices() // TODO ask Ajay
                        )
        {}

private:
    static std::vector< const Function<1>* > createFunctions(angle_t amplitude, angle_t phase) {
        std::vector< const Function<1>* > functions;
        functions.push_back(new ZeroFunction ); // x rotation
        functions.push_back(new ZeroFunction ); // y rotation

        // Z-rotation is the only degree of freedom for this mobilizer
        functions.push_back(new SinusoidFunction(amplitude, phase) ); // z rotation

        functions.push_back(new ZeroFunction ); // x translation
        functions.push_back(new ZeroFunction ); // y translation
        functions.push_back(new ZeroFunction ); // z translation

        return functions;        
    }

    static std::vector< std::vector<int> > createCoordIndices() {
        std::vector<int> oneCoordinateVec;
        oneCoordinateVec.push_back(0); // first and only generalized coordinate index

        std::vector< std::vector<int> > coordIndices;
        // All six functions take the one generalized coordinate
        coordIndices.push_back(std::vector<int>()); // 1, empty
        coordIndices.push_back(std::vector<int>()); // 2, empty
        coordIndices.push_back(oneCoordinateVec);   // 3, one coordinate, zero
        coordIndices.push_back(std::vector<int>()); // 4, empty
        coordIndices.push_back(std::vector<int>()); // 5, empty
        coordIndices.push_back(std::vector<int>()); // 6, empty

        return coordIndices;
    }
};



void testRiboseMobilizer() 
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);

    matter.setShowDefaultGeometry(false);

    // Put some hastily chosen mass there (doesn't help)
    Body::Rigid rigidBody;
    rigidBody.setDefaultRigidBodyMassProperties(MassProperties(
        mass_t(20.0*daltons),
        location_t(Vec3(0,0,0)*nanometers),
        Inertia(20.0)
        ));

    // One body anchored at C4 atom, 
    MobilizedBody::Weld c4Body( 
        matter.updGround(), 
        Rotation(-120*degrees, XAxis),
        rigidBody,
        Transform());
    // sphere for C4 atom
    decorations.addBodyFixedDecoration(
        c4Body.getMobilizedBodyIndex(), 
        Transform(),
        DecorativeSphere( length_t(0.5*angstroms) )
    );
    // sphere for C5 atom
    decorations.addBodyFixedDecoration(
        c4Body.getMobilizedBodyIndex(), 
        location_t(Vec3(-1.0,-1.0,0.5)*angstroms),
        DecorativeSphere( length_t(0.5*angstroms) )
    );

    decorations.addRubberBandLine(
        c4Body.getMobilizedBodyIndex(),
        Vec3(0),
        c4Body.getMobilizedBodyIndex(),
        location_t(Vec3(-1.0,-1.0,0.5)*angstroms),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));

    // One body anchored at C3 atom -- works
    // Pin version
    //MobilizedBody::Pin c3Body( 
    //    c4Body, 
    //    Transform(),
    //    rigidBody,
    //    Transform(location_t(Vec3(0,0,1.5)*angstroms))
    //    );

    // Function based pin version -- works
    //TestPinMobilizer c3Body( 
    //    c4Body, 
    //    Transform(),
    //    rigidBody,
    //    Transform(location_t(Vec3(0,0,1.5)*angstroms))
    //    );
    
    PseudorotationMobilizer c3Body( 
        c4Body, 
        Transform(),
        rigidBody,
        Transform(location_t(Vec3(0,0,1.5)*angstroms)),
        angle_t(36.4*degrees), // amplitude
        angle_t(-161.8*degrees) // phase
        );
    // sphere for C3 atom
    decorations.addBodyFixedDecoration(
        c3Body.getMobilizedBodyIndex(), 
        Transform(),
        DecorativeSphere( length_t(0.5*angstroms) )
    );
    // sphere for O3 atom
    decorations.addBodyFixedDecoration(
        c3Body.getMobilizedBodyIndex(), 
        location_t(Vec3(-1.0,1.0,-0.5)*angstroms),
        DecorativeSphere( length_t(0.5*angstroms) ).setColor(Vec3(1,0,0))
    );

    decorations.addRubberBandLine(
        c3Body.getMobilizedBodyIndex(),
        Vec3(0),
        c3Body.getMobilizedBodyIndex(),
        location_t(Vec3(-1.0,1.0,-0.5)*angstroms),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));
    decorations.addRubberBandLine(
        c4Body.getMobilizedBodyIndex(),
        Vec3(0),
        c3Body.getMobilizedBodyIndex(),
        Vec3(0),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));

    PseudorotationMobilizer c2Body( 
        c3Body, 
        Rotation( angle_t(-80*degrees), YAxis ),
        rigidBody,
        Transform(location_t(Vec3(0,0,1.5)*angstroms)),
        angle_t(35.8*degrees), // amplitude
        angle_t(-91.3*degrees) // phase
        );
    // sphere for C2 atom
    decorations.addBodyFixedDecoration(
        c2Body.getMobilizedBodyIndex(), 
        Transform(),
        DecorativeSphere( length_t(0.5*angstroms) )
    );
    // sphere for O2 atom
    decorations.addBodyFixedDecoration(
        c2Body.getMobilizedBodyIndex(), 
        location_t(Vec3(-1.0,1.0,-0.5)*angstroms),
        DecorativeSphere( length_t(0.5*angstroms) ).setColor(Vec3(1,0,0))
    );

    decorations.addRubberBandLine(
        c2Body.getMobilizedBodyIndex(),
        Vec3(0),
        c2Body.getMobilizedBodyIndex(),
        location_t(Vec3(-1.0,1.0,-0.5)*angstroms),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));
    decorations.addRubberBandLine(
        c3Body.getMobilizedBodyIndex(),
        Vec3(0),
        c2Body.getMobilizedBodyIndex(),
        Vec3(0),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));

    PseudorotationMobilizer c1Body( 
        c2Body, 
        Rotation( angle_t(-80*degrees), YAxis ),
        rigidBody,
        Transform(location_t(Vec3(0,0,1.5)*angstroms)),
        angle_t(37.6*degrees), // amplitude
        angle_t(52.8*degrees) // phase
        );
    // sphere for C1 atom
    decorations.addBodyFixedDecoration(
        c1Body.getMobilizedBodyIndex(), 
        Transform(),
        DecorativeSphere( length_t(0.5*angstroms) )
    );
    // sphere for N1 atom
    decorations.addBodyFixedDecoration(
        c1Body.getMobilizedBodyIndex(), 
        location_t(Vec3(-1.0,-1.0,-0.5)*angstroms),
        DecorativeSphere( length_t(0.5*angstroms) ).setColor(Vec3(0,0,1))
    );
    // sphere for O4 atom
    decorations.addBodyFixedDecoration(
        c1Body.getMobilizedBodyIndex(), 
        location_t(Vec3(1.0,0,-0.5)*angstroms),
        DecorativeSphere( length_t(0.5*angstroms) ).setColor(Vec3(1,0,0))
    );

    decorations.addRubberBandLine(
        c2Body.getMobilizedBodyIndex(),
        Vec3(0),
        c1Body.getMobilizedBodyIndex(),
        Vec3(0),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));
    decorations.addRubberBandLine(
        c1Body.getMobilizedBodyIndex(),
        Vec3(0),
        c1Body.getMobilizedBodyIndex(),
        location_t(Vec3(1.0,0,-0.5)*angstroms),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));
    decorations.addRubberBandLine(
        c1Body.getMobilizedBodyIndex(),
        Vec3(0),
        c1Body.getMobilizedBodyIndex(),
        location_t(Vec3(-1.0,-1.0,-0.5)*angstroms),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));
    decorations.addRubberBandLine(
        c4Body.getMobilizedBodyIndex(),
        Vec3(0),
        c1Body.getMobilizedBodyIndex(),
        location_t(Vec3(1.0,0,-0.5)*angstroms),
        DecorativeLine().setColor(Vec3(0,0,0)).setLineThickness(6));

    // Prescribed motion
    Constraint::ConstantSpeed(c3Body, 0.5);

    // Two constraint way works; one constraint way does not
    bool useTwoConstraints = false;

    if (useTwoConstraints) {
        // Constraints to make three generalized coordinates identical
        std::vector<MobilizedBodyIndex> c32bodies(2);
        c32bodies[0] = c3Body.getMobilizedBodyIndex();
        c32bodies[1] = c2Body.getMobilizedBodyIndex();
        std::vector<MobilizerQIndex> coordinates(2, MobilizerQIndex(0));
        Constraint::CoordinateCoupler(matter, new DifferenceFunction, c32bodies, coordinates);

        std::vector<MobilizedBodyIndex> c21bodies(2);
        c21bodies[0] = c2Body.getMobilizedBodyIndex();
        c21bodies[1] = c1Body.getMobilizedBodyIndex();
        Constraint::CoordinateCoupler(matter, new DifferenceFunction, c21bodies, coordinates);
    }
    else { // trying to get single constraint way to work
        // Try one constraint for all three mobilizers
        std::vector<MobilizedBodyIndex> c123Bodies(3);
        c123Bodies[0] = c1Body.getMobilizedBodyIndex();
        c123Bodies[1] = c2Body.getMobilizedBodyIndex();
        c123Bodies[2] = c3Body.getMobilizedBodyIndex();
        std::vector<MobilizerQIndex> coords3(3, MobilizerQIndex(0));
        Constraint::CoordinateCoupler(matter, new ThreeDifferencesFunction, c123Bodies, coords3);
    }

    system.updDefaultSubsystem().addEventReporter( new VTKEventReporter(system, 0.10) );

    system.realizeTopology();
    State& state = system.updDefaultState();
    
    // Simulate it.
    VerletIntegrator integ(system);
    //RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(50.0);
}

int main() 
{
    testRiboseMobilizer();

    return 0;
}
