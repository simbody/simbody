#ifndef __internalDynamics_hh__
#define __internalDynamics_hh__

//
// namespace with common defs
//

namespace InternalDynamics {

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

class Exception { 
public: 
    const char* mess; 
    Exception(const char* mess) : mess(mess) {}
};

};

#endif /* __internalDynamics_hh__ */
