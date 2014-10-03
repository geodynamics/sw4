###############################################################################
# CMake module to search for cencalvm library
#
# On success, the macro sets the following variables:
# CENCALVM_FOUND       = if the library found
# CENCALVM_LIBRARIES   = full path to the library
# CENCALVM_INCLUDE_DIR = where to find the library headers 
#
# Copyright (c) 2014 Eric Heien <emheien@ucdavis.edu>
#
###############################################################################

# Try to use OSGeo4W installation
if($ENV{CENCALVM_DIR})
    set(CENCALVM_HOME "$ENV{CENCALVM_DIR}") 
endif()

find_path(CENCALVM_INCLUDE_DIR cencalvm/query/VMQuery.h
    PATHS ${CENCALVM_HOME}/include
    DOC "Path to cencalvm library include directory")

find_library(CENCALVM_LIBRARY
    NAMES cencalvm
    PATHS ${CENCALVM_HOME}/lib
    DOC "Path to cencalvm library file")

find_library(ETREE_LIBRARY
    NAMES etree
    PATHS ${CENCALVM_HOME}/lib
    DOC "Path to etree library file")

if(CENCALVM_LIBRARY AND ETREE_LIBRARY)
    set(CENCALVM_LIBRARIES ${CENCALVM_LIBRARY} ${ETREE_LIBRARY})
endif(CENCALVM_LIBRARY AND ETREE_LIBRARY)

# Handle the QUIETLY and REQUIRED arguments and set CENCALVM_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CENCALVM DEFAULT_MSG
    CENCALVM_LIBRARIES
    CENCALVM_INCLUDE_DIR)
