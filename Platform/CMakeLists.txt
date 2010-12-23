#---------------------------------------------------
# Simbody/Prerequisites
#
# Simbody has dependencies on several low-level packages
# for which the details differ on different platforms.
# Dependencies are:
#    - pthreads
#    - Lapack
# Pthreads is standard on Linux platforms, and close enough
# on Macs. On Windows we have to supply our own.
# We provide our own high-performance Lapack libraries
# on all platforms, but you can provide your own if you
# prefer by giving the library names. However, the 
# signatures must match our include file SimTKlapack.h.
#
#----------------------------------------------------

cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

# TODO

