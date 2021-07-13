## FIND PROJ4 
if (NOT EXISTS "${PROJ4_DIR}")
  message(FATAL_ERROR "Given SEARCH_DIR does not exist: ${PROJ4_DIR}")
endif()

if (NOT IS_DIRECTORY "${PROJ4_DIR}")
    message(FATAL_ERROR "Given SEARCH_DIR is not a directory: ${PROJ4_DIR}")
endif()

# Find includes directory
find_path( PROJ4_INCLUDE_DIR proj_api.h
           PATHS  ${PROJ4_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

# Find libraries
find_library( PROJ4_LIBRARY NAMES proj
              PATHS ${PROJ4_DIR}/lib64
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PROJ4_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( PROJ4  DEFAULT_MSG
                                  PROJ4_INCLUDE_DIR
                                  PROJ4_LIBRARY )

if(NOT PROJ4_FOUND)
    message(FATAL_ERROR "PROJ4_DIR is not a path to a valid PROJ4 install")
endif()

message(STATUS "PROJ4 Includes: ${PROJ4_INCLUDE_DIR}")
message(STATUS "PROJ4 Libraries: ${PROJ4_LIBRARY}")

blt_import_library(NAME       proj4
                   TREAT_INCLUDES_AS_SYSTEM ON
                   DEFINES    HAVE_PROJ4=1
                   INCLUDES   ${PROJ4_INCLUDE_DIR}
                   LIBRARIES  ${PROJ4_LIBRARY}
                   EXPORTABLE ON)



## FIND LAPACK

# first Check for LAPACK_DIR
if (NOT EXISTS "${LAPACK_DIR}")
    message(FATAL_ERROR "Given LAPACK_DIR does not exist: ${LAPACK_DIR}")
endif()

if (NOT IS_DIRECTORY "${LAPACK_DIR}")
    message(FATAL_ERROR "Given LAPACK_DIR is not a directory: ${LAPACK_DIR}")
endif()

# Find includes directory
find_path( LAPACK_INCLUDE_DIR lapacke.h
           PATHS  ${LAPACK_DIR}/include/
                  ${LAPACK_DIR}/include/lua
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

# Find libraries
if (ENABLE_HIP)
         message(STATUS "Enable HIP is defined")
	set(lalibname "sci_cray" CACHE STRING "" FORCE)
else()
         message(STATUS "Enable HIP is not defined")
	set(lalibname "lapack" CACHE STRING "" FORCE)
endif()
message(STATUS "LAPACK looking for ${lalibname}")
find_library( LAPACK_LIBRARY NAMES ${lalibname}
              PATHS ${LAPACK_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LAPACK_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LAPACK  DEFAULT_MSG
                                  LAPACK_INCLUDE_DIR
                                  LAPACK_LIBRARY )

if(NOT LAPACK_FOUND)
    message(FATAL_ERROR "LAPACK_DIR is not a path to a valid LAPACK install")
endif()

message(STATUS "LAPACK Includes: ${LAPACK_INCLUDE_DIR}")
message(STATUS "LAPACK Libraries: ${LAPACK_LIBRARY}")

blt_import_library(NAME       lapack
                   TREAT_INCLUDES_AS_SYSTEM ON
                   DEFINES    HAVE_LAPACK=1
                   INCLUDES   ${LAPACK_INCLUDE_DIR}
                   LIBRARIES  ${LAPACK_LIBRARY}
                   EXPORTABLE ON)


set (SW4_HDR_FILES AllDims.h    ForcingTwilight.h            MaterialProperty.h      SecondOrderSection.h  boundaryConditionTypes.h  hip_policies.h
    AnisotropicMaterial.h       GeographicProjection.h       MaterialRfile.h         SfileOutput.h         caliper.h                 impose_cartesian_bc.h
    AnisotropicMaterialBlock.h  GridGenerator.h              MaterialSfile.h         Source.h              cf_interface.h            nearlyEqual.h
    Byteswapper.h               GridGeneratorGaussianHill.h  MaterialVolimagefile.h  StatMachine.h         csstime.h                 omp_policies.h
    CheckPoint.h                GridGeneratorGeneral.h       Mspace.h                SuperGrid.h           cuda_foralls.h            optimizedcuda.h
    CurvilinearInterface2.h     GridPointSource.h            Parallel_IO.h           TestEcons.h           cuda_policies.h           policies.h
    DataPatches.h               Image.h                      Patch.h                 TestEnergy.h          cuda_policies12.h         readhdf5.h
    ESSI3D.h                    Image3D.h                    Policies.h              TestLamb.h            curvilinear4sgc.h         sachdf5.h
    ESSI3DHDF5.h                MaterialBlock.h              Polynomial.h            TestPointSource.h     curvilinear4sgcX3.h       sacsubc.h
    EW.h                        MaterialData.h               Qspline.h               TestRayleighWave.h    f_interface.h             sacutils.h
    EtreeFile.h                 MaterialIfile.h              RandomizedMaterial.h    TestTwilight.h        foralls.h                 startEnd.h
    F77_FUNC.h                  MaterialInvtest.h            Require.h               TimeDep.h             forallsdecl.h             time_functions.h
    Filter.h                    MaterialPfile.h              Sarray.h                TimeSeries.h          hip_foralls.h             version.h 
    sw4.h)

set (SW4_SRC_FILES main.C EW.C Sarray.C version.C parseInputFile.C ForcingTwilight.C 
    curvilinearGrid.C boundaryOp.f bndryOpNoGhost.f90  bcfort.f twilightfort.f rhs4th3fort.f
    parallelStuff.C Source.C MaterialProperty.C MaterialData.C material.C setupRun.C 
    solve.C solerr3.f Parallel_IO.C Image.C GridPointSource.C MaterialBlock.C testsrc.f
    TimeSeries.C sacsubc.C SuperGrid.C addsgd.f velsum.f rayleighfort.f energy4.f TestRayleighWave.C 
    MaterialPfile.C Filter.C Polynomial.C SecondOrderSection.C time_functions.C Qspline.C 
    lamb_exact_numquad.f twilightsgfort.f EtreeFile.C MaterialIfile.C GeographicProjection.C 
    rhs4curvilinear.f curvilinear4.f rhs4curvilinearsg.f curvilinear4sg.f gradients.f Image3D.C 
    MaterialVolimagefile.C MaterialRfile.C randomfield3d.f innerloop-ani-sgstr-vc.f bcfortanisg.f
    AnisotropicMaterialBlock.C checkanisomtrl.f computedtaniso.f sacutils.C ilanisocurv.f 
    anisomtrltocurvilinear.f bcfreesurfcurvani.f tw_ani_stiff.f90 tw_aniso_force.f tw_aniso_force_tt.f
    updatememvar.f90 addmemvarforcing2.C addsg4wind.f90 consintp.C scalar_prod.f90 oddIoddJinterp.C evenIoddJinterp.C 
    oddIevenJinterp.C evenIevenJinterp.C CheckPoint.C Mspace.C RandomizedMaterial.C AllDims.C Patch.C ESSI3D.C 
    MaterialSfile.C MaterialInvtest.C geodyn.C ESSI3DHDF5.C sachdf5.C readhdf5.C CurvilinearInterface2.C 
    TestEcons.C TestTwilight.C  TestPointSource.C curvilinear4sgwind.C GridGeneratorGeneral.C GridGeneratorGaussianHill.C 
    GridGenerator.C RHS43DEV.C curvilinear4sgcX1.C SfileOutput.C addsgdc.C bcfortc.C bcfortanisgc.C bcfreesurfcurvanic.C 
    boundaryOpc.C energy4c.C checkanisomtrlc.C computedtanisoc.C curvilinear4sgc.C gradientsc.C randomfield3dc.C 
    innerloop-ani-sgstr-vcc.C ilanisocurvc.C rhs4curvilinearc.C rhs4curvilinearsgc.C rhs4th3fortc.C solerr3c.C testsrcc.C 
    tw_aniso_forcec.C tw_aniso_force_ttc.C velsumc.C twilightfortc.C twilightsgfortc.C tw_ani_stiffc.C
    anisomtrltocurvilinearc.C scalar_prodc.C updatememvarc.C addsg4windc.C bndryOpNoGhostc.C rhs4th3windc2.C rhs4th3windc.C
    quadpack/dqawse.f quadpack/d1mach.f quadpack/dqags.f quadpack/dqaws.f
    quadpack/dqc25s.f quadpack/dqmomo.f quadpack/dqpsrt.f quadpack/dqk15w.f quadpack/xerror.f quadpack/dqags.f quadpack/dqagse.f
    quadpack/dqc25s.f quadpack/dqcheb.f quadpack/dqwgts.f quadpack/dqk21.f quadpack/dqelg.f)



#
# All the CUDA statements
#
if (ENABLE_CUDA)
  blt_add_executable(NAME sw4
    HEADERS ${SW4_HDR_FILES}
    SOURCES ${SW4_SRC_FILES}
    DEPENDS_ON cuda mpi)
  blt_add_target_compile_flags(TO sw4 FLAGS -DENABLE_CUDA=1)
  set_property(TARGET sw4 PROPERTY CUDA_ARCHITECTURES 70)
  
  if (ENABLE_SPILL_WARNS)
    blt_add_target_compile_flags(TO sw4 FLAGS -Xptxas -v,--warn-on-spills)

    
  endif()
  #set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda ")
  blt_add_target_link_flags(TO sw4 FLAGS -lcuda)
endif()
#
# All the HIP statements
#
if (ENABLE_HIP)
  blt_add_executable(NAME sw4
    HEADERS ${SW4_HDR_FILES}
    SOURCES ${SW4_SRC_FILES}
    DEPENDS_ON hip mpi)
  blt_add_target_compile_flags(TO sw4 FLAGS -DENABLE_HIP=1)
endif()

target_link_libraries(sw4 PUBLIC RAJA)
target_link_libraries(sw4 PUBLIC umpire)
target_link_libraries(sw4 PUBLIC lapack)
target_link_libraries(sw4 PUBLIC gomp nvidia-ml)
target_link_libraries(sw4 PUBLIC proj4)
if (ENABLE_HDF5)
target_link_libraries(sw4 PUBLIC HDF5::HDF5)
endif()
#target_link_libraries(sw4 PUBLIC cuda)
set_property(TARGET sw4 PROPERTY LINKER_LANGUAGE CXX)
# Can't get target_link_libraries to add -lcuda. Just gets ignored. End run around the problem using BLT

