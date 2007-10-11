#include "SimTKsimbody.h"
#include <iostream>

using namespace std;
using namespace SimTK;

int main() {
    cout << "Testing... ";

    DuMMForceFieldSubsystem dumm;
    assert(dumm.getNAtoms() == 0);

    cout << "PASSED" << endl;

    return 0;
}

