
/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
 * Contributors:
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "SimTKmath.h"
#include "OptimizerRep.h"

namespace SimTK {
    OptimizerRep::~OptimizerRep() {

        if( jacDiff ) delete jacDiff;
        if( gradDiff ) delete gradDiff;
        if( cf ) delete cf;
        if( of ) delete of;

     }
    void OptimizerRep::setConvergenceTolerance( const Real tolerance ){

       convergenceTolerance = tolerance;
       return;
    }
    void OptimizerRep::setMaxIterations( const int iter ){

       maxIterations = iter;
       return;
    }
    void OptimizerRep::setLimitedMemoryHistory( const int history ){

       limitedMemoryHistory = history;
       return;
    }

    void OptimizerRep::setDiagnosticsLevel( const int  level ){

       diagnosticsLevel = level;
       return;
    }

    template<class T> bool setAdvancedOptionHelper(std::map<std::string,T> &optionMap, const std::string &option, const T &value){
        optionMap[option] = value;
        return true;
    }

    bool OptimizerRep::setAdvancedStrOption( const std::string &option, const std::string &value ){
        return setAdvancedOptionHelper(advancedStrOptions, option, value);
    }
    bool OptimizerRep::setAdvancedRealOption( const std::string &option, const Real value ){
        return setAdvancedOptionHelper(advancedRealOptions, option, value);
    }
    bool OptimizerRep::setAdvancedIntOption( const std::string &option, const int value ){
        return setAdvancedOptionHelper(advancedIntOptions, option, value);
    }
    bool OptimizerRep::setAdvancedBoolOption( const std::string &option, const bool value ){
        return setAdvancedOptionHelper(advancedBoolOptions, option, value);
    }

    template<class T> bool getAdvancedOptionHelper(const std::map<std::string,T> &optionMap, const std::string &option, T &value){
        typename std::map<std::string,T>::const_iterator iter = optionMap.find(option);
        if(iter != optionMap.end()) {
            value = iter->second;
            return true;
        } else return false;
    }

    bool OptimizerRep::getAdvancedStrOption( const std::string &option, std::string &value ) const{
        return getAdvancedOptionHelper(advancedStrOptions, option, value);
    }
    bool OptimizerRep::getAdvancedRealOption( const std::string &option, Real &value ) const{
        return getAdvancedOptionHelper(advancedRealOptions, option, value);
    }
    bool OptimizerRep::getAdvancedIntOption( const std::string &option, int &value ) const{
        return getAdvancedOptionHelper(advancedIntOptions, option, value);
    }
    bool OptimizerRep::getAdvancedBoolOption( const std::string &option, bool &value ) const{
        return getAdvancedOptionHelper(advancedBoolOptions, option, value);
    }


    void OptimizerRep::useNumericalGradient( const bool flag ) {

       if( flag ) {     // turn on numerical gradients
           if( !numericalGradient ) initNumericalGrad();       
           numericalGradient = true;
       } else {        // turn off numerical graidents
           numericalGradient = false;
       }
       return;
    }
    void OptimizerRep::useNumericalJacobian( const bool flag ) {

       if( flag ) {     // turn on numerical jacbobian
           if( !numericalJacobian )  initNumericalJac();       
           numericalJacobian = true;
       } else {
           numericalJacobian = false;
       }
       return;
    }

    void OptimizerRep::initNumericalJac() {  // instaniates a jacobian Differentiator

        cf      = new SysConstraintFunc(sysp->numConstraints, sysp->numParameters, sysp );
        jacDiff = new Differentiator(*cf);  // construct Differentiator

    }
    void OptimizerRep::initNumericalGrad() {  // instaniates a gradient Differentiator

        of       = new SysObjectiveFunc( sysp->numParameters, sysp );
        gradDiff = new Differentiator(*of );  // construct Differentiator

     }


} // namespace SimTK
