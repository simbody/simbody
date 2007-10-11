#include "SimTKsimbody.h"
#include <iostream>

using namespace std;
using namespace SimTK;

int main() {
    cout << "Testing... ";

    DuMMForceFieldSubsystem dumm;

    DuMM::AtomId atomId = dumm.addAtom(DuMM::InvalidChargedAtomTypeId);

    cout << "PASSED" << endl;

    return 0;
}

