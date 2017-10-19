# This file was copied from Eigen.
# Eigen is primarily MPL2 licensed. See COPYING.MPL2 and these links:
#   http://www.mozilla.org/MPL/2.0/
#     http://www.mozilla.org/MPL/2.0/FAQ.html
if (ADOLC_INCLUDES AND ADOLC_LIBRARIES)
  set(ADOLC_FIND_QUIETLY TRUE)
endif (ADOLC_INCLUDES AND ADOLC_LIBRARIES)

find_path(ADOLC_INCLUDES
  NAMES
  adolc/adtl.h
  PATHS
  $ENV{ADOLCDIR}/include
)

find_library(ADOLC_LIBRARIES adolc PATHS $ENV{ADOLCDIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ADOLC DEFAULT_MSG
                                  ADOLC_INCLUDES ADOLC_LIBRARIES)

mark_as_advanced(ADOLC_INCLUDES ADOLC_LIBRARIES)
