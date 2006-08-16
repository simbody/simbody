#include "Simmath.h"
#include "optimizer.h"
#include "optimizerImplementation.h"

 namespace SimTK {
smOptimizer::smOptimizer(int dimension) {
        optimizerImplementation *optPtr =  new optimizerImplementation(dimension);
        data = (void *)optPtr;

}
} // namespace SimTK
