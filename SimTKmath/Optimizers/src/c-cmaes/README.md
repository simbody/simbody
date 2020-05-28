c-cmaes
=======

CMA-ES written in ANSI C in a fairly object-oriented style. 

For the general purpose of this software see `doc.txt` or 
[here](https://www.lri.fr/~hansen/cmaesintro.html), 
for more documentation on this library see `docfunctions.txt`, 
for how to start see below.

------------------------------------------------------------------

Files in this Repository
------------------------

- README.md : this file
- LICENSE : users agreement (no worries)
- doc.txt : describes general purpose and an application issue
- docfunctions.txt : Documentation of the library functions. 
- cmaes_interface.h : User interface header.   
- example_short.c : Very short example source code. The purpose of
               the example codes is to be edited/extended.  
- example_restarts.c : implements additional restarts with increasing
               population size (Auger & Hansen 2005). 
- example_boundary.c : combination with boundary handling (box 
     constraints)
- example_noise.c : (future versions) implements an additional 
               uncertainty handling (Hansen et al 2009). 
- cmaes.h : Header, e.g. declaration of struct cmaes_t.  
- cmaes.c : Source code.
- cmaes_initials.par : Parameters to be read by the cmaes, e.g. problem
      dimension, initial search point and initial standard deviation. 
      This file should be edited. 
- cmaes_signals.par : File for controlling the running program. Printing 
      or writing to a file/console can be set on/off while the program 
      is running. Regular termination can be forced. On delivery
      the writing is in accordance with the plotting using: 
- plotcmaesdat.m : Plots default output files in Matlab or Octave.
- plotcmaesdat.sci : Plots default output files in Scilab. 
- boundary_transformation.c : implements a boundary transformation
- boundary_transformation.h : header file. 

Files You May Need to Edit
--------------------------

- example_*.c:  Plug in the objective function (pointer) that should 
    be minimized. 
- cmaes_initials.par: Parameter file for changing e.g. initial values and
    stopping criteria without recompiling. 
- cmaes_signals.par: File to control termination and output during 
    runtime. 


Output files written by cmaes_t
-------------------------------

- actparcmaes.par : Parameters as actually used by the program. The
    actual parameter setting is appended to the file after each start 
    of the cmaes. 
- errcmaes.err : Error messages. 


------------------------------------------------------------------

HOW TO START
============

  0) get code via `git ...` or download button
  
  A1) Take five minutes to look at file `example_short.c`. 

  A2) You might have a glance at the documentation provided in file
     `docfunctions.txt`.

  A3) You might have a glance at `cmaes_initials.par`, where input parameters
     are defined. 

  B1) Compile and run the example program. Compilation e.g. with 
     the GNU c-compiler in the `src` folder:

	gcc -Wall -o evo cmaes.c example_short.c -lm
     
  and run with `evo` or `./evo`. Take a look at the output. 

  B2a) (optional but highly recommended: plotting) Invoke Scilab (freely 
     available for Linux/Windows/Mac) or Matlab/Octave, change to the 
     working directory and type (Scilab)
        `getf('plotcmaesdat.sci'); plotcmaesdat;` 
     or (Matlab/Octave)
        `plotcmaesdat;`
     You need to have the file `plotcmaesdat.sci` or `.m` and the
     output data files in the working directory. You get a nice plot
     of the executed run.
     The same works with [`cma.py`](https://pypi.python.org/pypi/cma) via
       `python cma.py plot`

  B2b) (optional) Change (increase) problem dimension and/or problem
     number in file initials.par and re-run.

  B2c) (optional) Change problem dimension in initials.par to 300 and
     change output verbosity via file signals.par while the program
     is running: change e.g. "print fewinfo 200" into "print fewinfo
     -200" *and back*. Read comments. 

  B2d) Change back problem dimension.  

  5) Now you are ready to inspect and edit `example_restarts.c` or `example_boundary.c`
    to plug in the function you want to optimize. Refer to `doc.txt` and [see here](https://www.lri.fr/~hansen/cmaes_inmatlab.html#practical) for
    a practical issue on objective function design. Refer to
    `docfunctions.txt` to find more documentation about the functions in 
    this package. 

  6) Check "obligatory settings" part in `initials.par` regarding your
     function. Make sure that the scale of all objective parameter
     components of the function is somewhat similar and sigma
     corresponds to about 1/4 of the respective search intervals.

  7) output files are overwritten with each program call.  


Questions? [go here](https://github.com/cma-es/c-cma-es/issues/new) or 
send an email to hansen at lri dot fr. 

See also: 
- [Practical hints](https://www.lri.fr/~hansen/cmaes_inmatlab.html#practical)
- [CMA-ES short intro](https://www.lri.fr/~hansen/cmaesintro.html)
- [CMA-ES on wikipedia](http://en.wikipedia.org/wiki/CMA-ES)

