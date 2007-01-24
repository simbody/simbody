#include <iostream>
#include "SimTKSimbody.h"
#include "chemistry/Atom.h"
#include "AtomModel.h"
#include "MoleculeModeler.h"
#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/VTKReporter.h"

using namespace std;
using namespace SimTK;

// This is a test program to guide implementation of the molecular
// modeling layer that has yet to be constructed.
int main()
{
	Vec3 v1(-0.5, 0.1, 0);
	cout << v1 << endl;

	MoleculeModeler modeler;
	
	Atom argonAtom(ChemicalElement::Argon);
	cout << argonAtom.getDefaultPosition() << endl;
	argonAtom.setDefaultPosition(Vec3(9,9,8));
	cout << argonAtom.getDefaultPosition() << endl;

	const AtomModel argonModel1 = modeler.addAtomLike(argonAtom, Vec3(-0.5, 0.1, 0));
	const AtomModel argonModel2 = modeler.addAtomLike(argonAtom, Vec3( 0.5, 0, 0));
	const AtomModel argonModel3 = modeler.addAtomLike(argonAtom, Vec3( 0.0, -0.5, 0.1));
	const AtomModel argonModel4 = modeler.addAtomLike(argonAtom, Vec3( 0.0, -0.1, 0.7));

	SimTK::Vec3 vec1 = argonModel1.getDefaultPosition();
	cout << "  Position 1 = " << vec1 << endl;
	cout << "  Position 2 = " << argonModel2.getDefaultPosition() << endl;
	cout << "  Position 3 = " << argonModel3.getDefaultPosition() << endl;
	cout << "  Position 4 = " << argonModel4.getDefaultPosition() << endl;	

	State state;
	const MultibodySystem & system = modeler.realizeSystem(state);

	VTKReporter vtk(system);

	RungeKuttaMerson study(system, state);
	Real stepInterval = 0.05; // picoseconds
	Real time = 0;
	state.updTime() = time;

	study.initialize();

	while (time < 1000 * stepInterval)
	{
		study.step(time);

		cout << "Energy = " << system.getEnergy(state) << endl;
		cout << "  Position 1 = " << argonModel1.getPosition(state) << endl;
		cout << "  Position 2 = " << argonModel2.getPosition(state) << endl;
		cout << endl; // blank line

		vtk.report(state);

		time += stepInterval;
	}

	return 0;	
}
