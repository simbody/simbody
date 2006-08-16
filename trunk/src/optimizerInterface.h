#ifndef _SimTK_OPTIMIZER_INTERFACE_H_
#define _SimTK_OPTIMIZER_INTERFACE_H_

 namespace SimTK {

class smOptimizerInterface { 

public:
    virtual ~smOptimizerInterface() {};
   

    /* by default class chooses algorithm based on size of problem. 
       The algorithm can be overridden by setting the Algo param. */
       /* should there be a set method for each parameter ? */
    virtual int setOptimizerParameters(int, double *) = 0;
    virtual int getOptimizerParameters(int, double *) = 0;
    virtual int setObjectiveFunction( int *) = 0; // TODO use costFunction
    virtual int setInitialPosition( double * ) = 0; // TODO use simmatrix
    virtual int optimize(double *) = 0; // checks to see if space needs to be allocated
                             // constructor sets internal flag to require allocation
                             // if dimension or algorithm changes the reallocate flag is
                             // set to trigger freeing the old mem and allocating new.


}; // end class optimizeInterface
} // namespace SimTK
#endif  //_SimTK_OPTIMIZER_INTERFACE_H_
