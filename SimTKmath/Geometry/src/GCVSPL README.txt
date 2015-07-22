C version of GCVSPL

This is a version of Herman Woltring's classic GCVSPL fortran program
that has been translated to c. Herman's original fortran code has been
translated with the free fortran to c translator, f2c. (For more
information on f2c consult netlib at research.att.com where f2c is
archived.)

I have set f2c to translate the program to ANSI c. The test program
included here is a hand coded version of Herman's original test program.
This is included to show how gcvspl should be accessed in c. Note that
the code produced by f2c is almost unreadable but it has compiled and
run for me on a number of different operating systems (Unix,DOS,Mac).
I did the hand coded example because even though the f2c conversion of
the original worked it was unreadable.

When you use gcvspl you with need the f2c.h file included here. Only a small
number of the data types and structures defined in this file are actually
used in the gcvspl.c code, so you can remove alot of the content of this
file if you like.

Lastly, there is an interesting addition to gcvspl that Ton van den
Bogert mentioned awhile back. You can set the cutoff frequency of gcvspl
or inquire what it selected as the cutoff frequency when it is run in
automatic mode. This is not demonstrated in this code. However, I
have implemented both these options in my gait package in the
interface code to gcvspl (in the file filter.for of the program ANZ).
Look in the same directory as this archive for more info on how get
those files.

dwight meglan
