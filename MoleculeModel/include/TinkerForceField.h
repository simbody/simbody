#ifndef TINKERFORCEFIELD_H
#define TINKERFORCEFIELD_H

#include "SimTKSimbody.h"

enum VdwType {LENNARD-JONES};
enum RadiusRule {ARITHMETIC};
enum RadiusType {R-MIN};
enum RadiusSize {RADIUS};
enum EpsilonRule {GEOMETRIC};

class TinkerForceField {
private:
	SimTK::String name;
	VdwType vdwType;
	RadiusRule radiusRule;
	RadiusType radiusType;
	RadiusSize radiusSize;
	EpsilonRule epsilonRule;
	double vdw-14-scale;
	double chg-14-scale;
	double dielectric;
};

#endif /*TINKERFORCEFIELD_H*/
