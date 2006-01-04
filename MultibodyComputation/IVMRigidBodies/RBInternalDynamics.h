#ifndef __internalDynamics_hh__
#define __internalDynamics_hh__

#include "cdsList.h"
#include "cdsString.h"

//
// namespace with common defs
//

namespace IVMInternalDynamics {

enum {                            //possible values for verbose (bit field)
    printCoords          = 1,
    printResetCM         = 2,
    printVelFromCartCost = 4,
    printTemperature     = 8,
    printEnergy          = 16,
    printCMVel           = 32,
    printNodeForce       = 64,
    printNodePos         = 128,
    printNodeTheta       = 256,
    printStepDebug       = 512,   //integrator/minimizer debug statements
    printStepInfo        = 1024,  //integrator/minimizer info  statements
    printNodeDef         = 2048,
    printLoopDebug       = 4096,
    printLoopInfo        = 8192
};

struct IVMHingeSpec {
    CDSString type;
    int atom0;  //extra data
    int atom1;
    IVMHingeSpec(const char* type="none") : type(type), atom0(-2), atom1(-2) {}
    CDSList<int> aList;
};

//  ostream& operator<<(      ostream&  s,
//            const HingeSpec& h);

class IVMException { 
public: 
    const char* mess; 
    IVMException(const char* mess) : mess(mess) {}
};
   
};

#endif /* __internalDynamics_hh__ */
