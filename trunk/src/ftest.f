      program main
      external  fsmmallocoptimizer
      external  fsmsetoptparams
      external  fsmsetcostfunc
      external  fsmdumpoptimizerstate
      external  fsmstartopt
      external  fsmgetoptimizerresult
      external  fsmfreeopt
      external  costfunction
      integer*8 handle
      integer   status
      double precision val
      double precision result(2),initialvalue(2)
      

      call fsmmallocoptimizer("L", 2, handle, status) 
      if( status .NE. 0 ) write(*,*)"malloc failed status=",status

      val = 100
      call fsmsetoptparams(handle,"FUNCTION EVALUATIONS", val, status) 
      if( status .NE. 0 ) write(*,*)"number func eval failed status=",
     * status

      val = 1.0
      call fsmsetoptparams(handle,"STEP LENGTH", val, status) 
      if( status .NE. 0 ) write(*,*)"step length failed status=",status

      val = 0.0001; 
      call fsmsetoptparams(handle,"GRADIENT", val, status) 
      if( status .NE. 0 ) write(*,*)"gradient failed status=",status

      val =0.9; 
      call fsmsetoptparams(handle,"ACCURACY",val, status) 
      if( status .NE. 0 ) write(*,*)"accuracy failed status=",status

      initialvalue(1) =  100 
      initialvalue(2) = -100 
      call fsmsetoptparams( handle, "INITAL VALUES", initialvalue,
     * status )
      if( status .NE. 0 ) write(*,*)"initial values fail status=",status

      call fsmsetcostfunc( handle, costfunction, status )
      if( status .NE. 0 ) write(*,*)"cost func fail status=",status

      call fsmdumpoptimizerstate(handle, status) 
      if( status .NE. 0 ) write(*,*)"dump state fail status=",status

      call fsmstartopt(handle, status) 
      if( status .NE. 0 ) write(*,*)"start opt failed status=",status

      call fsmgetoptimizerresult(handle,result,status) 
      if( status .NE. 0 ) write(*,*)"get results failed status=",status

      write(*,*)"RESULTS=",result(1),result(2)

      call fsmfreeopt(handle, status );
      end

      subroutine costfunction( pos, f, g)
      double precision pos(2), f, g(2)

      double precision x, y

      x = pos(1)
      y = pos(2)
   
      f    = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y
      g(1) = 3*x + 2*y -2
      g(2) = 2*x + 6*y +8

c      write(*,*)"F costfunction pos=",pos(1),pos(2)," f=",f,
c     * "  g=",g(1)," ",g(2)

      end
