if (ENABLE_CUDA)
  blt_add_executable(NAME sw4
    HEADERS AllDims.h                   ForcingTwilight.h            MaterialProperty.h      SecondOrderSection.h  boundaryConditionTypes.h  hip_policies.h
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
    sw4.h
    SOURCES EW.C Sarray.C version.C parseInputFile.C ForcingTwilight.C 
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
    GridGenerator.C RHS43DEV.C curvilinear4sgcX1.C SfileOutput.C
    DEPENDS_ON cuda mpi)
  blt_add_target_compile_flags(TO sw4 FLAGS -DENABLE_CUDA=1)
  set_property(TARGET sw4 PROPERTY CUDA_ARCHITECTURES 70)
  target_link_libraries(sw4 PUBLIC RAJA)
  if (ENABLE_SPILL_WARNS)
    blt_add_target_compile_flags(TO sw4 FLAGS -Xptxas -v,--warn-on-spills)
  endif()
  #set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda ")
endif()
