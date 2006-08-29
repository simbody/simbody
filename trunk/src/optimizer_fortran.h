
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
#undef F77NAME
#undef F77FUNC

#define FSM_MALLOCOPTIMIZER     F77NAME(FSMMALLOCOPTIMIZER,    fsmmallocoptimizer)
#define FSM_FREEOPT             F77NAME(FSMFREEOPTIMIZER,            fsmfreeoptimizer)
#define FSM_SETOPTPARMS         F77NAME(FSMSETOPTIMIZERPARAMETERS,        fsmsetoptimizerparameters)
#define FSM_SETCOSTFUNC         F77NAME(FSMSETCOSTFUNCTION,        fsmsetcostfunction)
#define FSM_DUMPOPTIMIZERSTATE  F77NAME(FSMDUMPOPTIMIZERSTATE, fsmdumpoptimizerstate)
#define FSM_RUNOPT              F77NAME(FSMRUNOPTIMIZER,           fsmrunoptimizer)

#if defined(SIMTK_FORTRAN_UPPERCASE)

#define F77NAME( uname, lname)     F77FUNC(uname)

#else

#define F77NAME( uname, lname)     F77FUNC(lname)

#endif

#if defined(SIMTK_FORTRAN_ONE_UNDERSCORE)

#define F77FUNC(name) name##_

#elif defined(SIMTK_FORTRAN_TWO_UNDERSCORE) 

#define F77FUNC(name) name##__

#else

#define F77FUNC(name) name

#endif

typedef long* FORTRAN_HANDLE;

extern "C" void FSM_FREEOPT(FORTRAN_HANDLE);
extern "C" void FSM_MALLOCOPTIMIZER( int*, FORTRAN_HANDLE,  smStatus* );
extern "C" void FSM_DUMPOPTIMIZERSTATE( FORTRAN_HANDLE, smStatus* );
extern "C" void FSM_SETOPTPARMS( FORTRAN_HANDLE, char *, double*, smStatus* );
extern "C" void FSM_SETCOSTFUNC(FORTRAN_HANDLE, void(costFunction)(int, double *, double*, double*, void *), smStatus*);
extern "C" void FSM_RUNOPT( FORTRAN_HANDLE handle, double *, smStatus*);



void FSM_SETCOSTFUNC( FORTRAN_HANDLE handle, void(costFunction)(int, double *, double*, double*, void *), smStatus *status){

    *status = ((SimTK::OptimizerImplementation *)((long)*handle))->setObjectiveFunction(costFunction);
}

void FSM_SETOPTPARMS( FORTRAN_HANDLE handle, char *parameter, double *values, smStatus *status){

   unsigned int param;

   param = ((SimTK::OptimizerImplementation *)((long)*handle))->optParamStringToValue( parameter );

   *status = ((SimTK::OptimizerImplementation *)((long)*handle))->setOptimizerParameters(param,values);
   return;
}
void FSM_RUNOPT( FORTRAN_HANDLE handle, double *results, smStatus* status){

    *status = ((SimTK::OptimizerImplementation *)((long)*handle))->optimize(results);
    return;
}
void FSM_MALLOCOPTIMIZER( int *n, FORTRAN_HANDLE handle, smStatus *status){

      *status = SUCCESS; 
      SimTK::OptimizerImplementation *opt = new SimTK::OptimizerImplementation(*n);
      *handle = (long)opt;
      return;
}

void FSM_FREEOPT( FORTRAN_HANDLE handle){

    delete ((SimTK::OptimizerImplementation *)((long)*handle));
    return;
}
void FSM_DUMPOPTIMIZERSTATE(FORTRAN_HANDLE  handle, smStatus *status) {

    *status = smDumpOptimizerState( (void *)((long)*handle) );
    return;
}

