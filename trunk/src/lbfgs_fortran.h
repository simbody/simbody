
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

extern "C" void FSM_FREEOPT(FORTRAN_HANDLE);
extern "C" void FSM_MALLOCOPTIMIZER( int*, FORTRAN_HANDLE,  smStatus* );
extern "C" void FSM_DUMPOPTIMIZERSTATE( FORTRAN_HANDLE, smStatus* );
extern "C" void FSM_SETOPTPARMS( FORTRAN_HANDLE, char *, double*, smStatus* );
extern "C" void FSM_SETCOSTFUNC(FORTRAN_HANDLE, void(costFunction)(double *, double*, double*), smStatus*);
extern "C" void FSM_RUNOPT( FORTRAN_HANDLE handle, double *, smStatus*);



void FSM_SETCOSTFUNC( FORTRAN_HANDLE handle, void(costFunction)(double *, double*, double*), smStatus *status){

    *status =  smSetCostFunction( (void *)((long)*handle), costFunction );
}

void FSM_SETOPTPARMS( FORTRAN_HANDLE handle, char *parameter, double *values, smStatus *status){

   unsigned int param;

   if( 0 == strncmp( "FUNCION_EVALUATIONS", parameter, 1) ) {
     param = MAX_FUNCTION_EVALUATIONS;
   } else if( 0 == strncmp( "STEP_LENGTH", parameter, 1)) {
     param = DEFAULT_STEP_LENGTH;
   } else if( 0 == strncmp( "INITIAL_VALUES", parameter, 1)) {
     param = INITIAL_VALUES;
   } else if( 0 == strncmp( "TOLERANCE", parameter, 1)) {
     param = TRACE;
   } else if( 0 == strncmp( "GRADIENT", parameter, 1)) {
     param = GRADIENT_CONVERGENCE_TOLERANCE;
   } else if( 0 == strncmp( "ACCURACY", parameter, 1)) {
     param = LINE_SEARCH_ACCURACY;
   }

   *status =  smSetOptimizerParameters( (void *)((long)*handle), param, values );
   return;
}
void FSM_RUNOPT( FORTRAN_HANDLE handle, double *results, smStatus* status){

    *status = smRunOptimizer( (void *)((long)*handle), results );
    return;
}
void FSM_MALLOCOPTIMIZER( int *n, FORTRAN_HANDLE handle, smStatus *status){

      int opt_type;

      *handle = (long)smMallocOptimizer( *n, status);
      return;
}

void FSM_FREEOPT( FORTRAN_HANDLE handle){

    smFreeOptimizer( (void *)((long)*handle) );
    return;
}
void FSM_DUMPOPTIMIZERSTATE(FORTRAN_HANDLE  handle, smStatus *status) {

    *status = smDumpOptimizerState( (void *)((long)*handle) );
    return;
}

