#include <iostream>
#include "SimTKSimbody.h"
#include "chemistry/Atom.h"
#include "AtomModel.h"
#include "MoleculeModeler.h"
#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/VTKReporter.h"

using namespace SimTK;

int main() {
	MoleculeModeler modeler;

	State state;
	const MultibodySystem & system = modeler.realizeSystem(state);

	return 0;
}
