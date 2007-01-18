#ifndef TINKERATOMCLASS_H_
#define TINKERATOMCLASS_H_

class TinkerAtomClass {
private:
	SimTK::String name;
	int index;
	AmberAtomType amberAtomType;
	int nBonds;
	ChemicalElement element;
	
	double vanDerWaalsRadius;
	double vanDerWaalsWellDepth;
	
	double charge;
	
	// double mass; // If we want to precisely reproduce the forcefield
};

#endif /*TINKERATOMCLASS_H_*/
