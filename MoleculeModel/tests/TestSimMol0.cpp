#include <iostream>
#include "Simbody.h"
#include "OxygenMolecule.h"
#include "MoleculeModeler.h"

using namespace std;
using namespace SimTK;

// This is a test program to guide implementation of the molecular
// modeling layer that has yet to be constructed.
int main() 
{
	MoleculeModeler modeler;
	OxygenMolecule oxygen;
	modeler.addMolecule(oxygen);
	MultibodySystem system = modeler.system();
	
	State state;
	RungeKuttaMerson study(system, state);
	
	Real stepInterval = 0.05; // picoseconds
	Real time = stepInterval;
	while (time < 20 * stepInterval) 
	{
		study.step(time);
		
		cout << "Energy = " << system.getEnergy(state) << endl;
		
		time += stepInterval;
	}

	return 0;	
}
