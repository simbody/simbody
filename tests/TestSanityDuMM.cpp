#include "SimTKsimbody.h"
#include <iostream>

using namespace std;
using namespace SimTK;

int main() {
try {

    cout << "Testing... ";

    DuMMForceFieldSubsystem dumm;
    assert(dumm.getNAtoms() == 0);

    dumm.dumpCForceFieldParameters(cout);

    cout << "PASSED" << endl;

    return 0;

}
catch (...) {
    cout << "FAILED" << endl;
    return 1;
}

}

