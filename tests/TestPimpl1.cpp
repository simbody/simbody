
#include "SimTKsimbody.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#ifdef _MSC_VER
    extern template class SimTK::PIMPLHandle<SimTK::MobilizedBody, SimTK::MobilizedBodyImpl>;
#endif

#include <iostream>

int main() 
{
    std::cout << "passed" << std::endl;

    return 0;
}
