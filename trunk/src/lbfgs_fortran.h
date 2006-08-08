
#undef F77NAME
#undef F77FUNC

#define FSM_MALLOCOPTIMIZER     F77NAME(FSMMALLOCOPTIMIZER,    fsmmallocoptimizer)
#define FSM_FREEOPT             F77NAME(FSMFREEOPT,            fsmfreeopt)
#define FSM_SETOPTPARMS         F77NAME(FSMSETOPTPARMS,        fsmsetoptparams)
#define FSM_SETCOSTFUNC         F77NAME(FSMSETCOSTFUNC,        fsmsetcostfunc)
#define FSM_DUMPOPTIMIZERSTATE  F77NAME(FSMDUMPOPTIMIZERSTATE, fsmdumpoptimizerstate)
#define FSM_GETOPTIMIZERRESULT  F77NAME(FSMGETOPTIMIZERRESULT, fsmgetoptimizerresult)
#define FSM_STARTOPT            F77NAME(FSMSTARTOPT,           fsmstartopt)

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

extern "C" void FSM_FREEOPT(FORTRAN_HANDLE, smStatus*);
extern "C" void FSM_MALLOCOPTIMIZER( char*, int*, FORTRAN_HANDLE,  smStatus* );
extern "C" void FSM_DUMPOPTIMIZERSTATE( FORTRAN_HANDLE, smStatus* );
extern "C" void FSM_SETOPTPARMS( FORTRAN_HANDLE, char *, double*, smStatus* );
extern "C" void FSM_SETCOSTFUNC(FORTRAN_HANDLE, void(costFunction)(double *, double*, double*), smStatus*);
extern "C" void FSM_GETOPTIMIZERRESULT( FORTRAN_HANDLE, double *, smStatus* );
extern "C" void FSM_STARTOPT( FORTRAN_HANDLE handle, smStatus*);


void  FSM_GETOPTIMIZERRESULT( FORTRAN_HANDLE handle, double *results, smStatus *status){

    *status = smGetOptimizerResults( (void *)((long)*handle), results );
    return;
}

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
void FSM_STARTOPT( FORTRAN_HANDLE handle, smStatus* status){

    *status = smStartOptimizer( (void *)((long)*handle) );
    return;
}
void FSM_MALLOCOPTIMIZER( char *type, int *n, FORTRAN_HANDLE handle, smStatus *status){

      int opt_type;

      if( 0 == strncmp( "LBFGS", type, 1 )) {
         opt_type = LBFGS;
      }
      
      *handle = (long)smMallocOptimizer(opt_type, *n, status);
      return;
}

void FSM_FREEOPT( FORTRAN_HANDLE handle, smStatus *status){

    *status = smFreeOptimizer( (void *)((long)*handle) );
    return;
}
void FSM_DUMPOPTIMIZERSTATE(FORTRAN_HANDLE  handle, smStatus *status) {

    *status = smDumpOptimizerState( (void *)((long)*handle) );
    return;
}

