      program main
      external  fsmMallocOptimizer
      external  fsmSetoptimizerParameters
      external  fsmSetCostFunction
      external  fsmDumpOptimizerState
      external  fsmRunOptimizer
      external  fsmFreeOptimizer
      external  costfunction
      integer*8 handle
      double precision val
      double precision result(2)
      

      call fsmMallocOptimizer( 2, handle) 

      val = 100
      call fsmSetOptimizerParameters(handle,"FUNCTION EVALUATIONS", val,
     * ) 

      val = 1.0

      call fsmSetOptimizerParameters(handle,"STEP LENGTH", val) 

      val = 0.0001; 
      call fsmSetOptimizerParameters(handle,"GRADIENT", val) 

      val =0.9; 
      call fsmSetOptimizerParameters(handle,"ACCURACY",val) 

      result(1) =  100 
      result(2) = -100 

      call fsmSetCostFunction( handle, costfunction )

      call fsmDumpOptimizerState(handle) 

      call fsmRunOptimizer(handle, result) 
      write(*,*)"RESULTS= ",result(1),result(2)

      call fsmFreeOptimizer(handle );
      end

      subroutine costfunction( n, pos, f, g, user_data )
      integer n
      double precision pos(2), f, g(2), user_data

      double precision x, y

      x = pos(1)
      y = pos(2)
   
      f    = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y
      g(1) = 3*x + 2*y -2
      g(2) = 2*x + 6*y +8

c      write(*,*)"F costfunction pos=",pos(1),pos(2)," f=",f,
c     * "  g=",g(1)," ",g(2)

      end
