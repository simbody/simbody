#include <iostream>
#include "Simbody.h"
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
	MoleculeModeler modeler;
	
	Atom argonAtom(ChemicalElement::Argon);
	AtomModel argonModel1 = modeler.addAtomLike(argonAtom, Vec3(-0.5, 0.1, 0));
	AtomModel argonModel2 = modeler.addAtomLike(argonAtom, Vec3( 0.5, 0, 0));
	modeler.addAtomLike(argonAtom, Vec3( 0.0, -0.5, 0.1));
	modeler.addAtomLike(argonAtom, Vec3( 0.0, -0.1, 0.7));
	
	State state;
	MultibodySystem & system = modeler.realizeSystem(state);

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
