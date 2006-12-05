#include <iostream>
#include "Simbody.h"
#include "Protein.h"
#include "Ethane.h"
#include "MoleculeModeler.h"

// TODO - include this in Simbody.h
#include "simbody/internal/NumericalMethods.h"

using namespace std;
using namespace SimTK;

// This is a test program to guide implementation of the molecular
// modeling layer that has yet to be constructed.
int main() 
{
	Protein protein("ADA");
	Ethane ethane;
	MoleculeModeler modeler;
	modeler.addMolecule(protein).addMolecule(ethane.setPosition(1,2,3));
	MultibodySystem system = modeler.system();
	
	State state;
	RungeKuttaMerson study(system, state);
	
	Real stepInterval = 0.05; // picoseconds
	Real time = stepInterval;
	while (time < 20 * stepInterval) 
	{
		study.step(time);
		
		cout << "Energy = " << system.getEnergy(state) << endl;
		
		cout 
			<< "Psi torsion of residue 1 = " 
			<< modeler.getPsiTorsion(state, protein.residue(1)) 
			<< endl;
			
		time += stepInterval;
	}

	return 0;	
}
