# Native and RAJA kernels retain their existing source manifests.
set(SW4_NATIVE_CXX_SOURCES src/EW.C src/Sarray.C src/version.C src/parseInputFile.C
    src/ForcingTwilight.C src/curvilinearGrid.C src/parallelStuff.C src/Source.C
    src/MaterialProperty.C src/MaterialData.C src/material.C src/setupRun.C
    src/solve.C src/Parallel_IO.C src/Image.C src/GridPointSource.C src/MaterialBlock.C
    src/TimeSeries.C src/sacsubc.C src/SuperGrid.C src/TestRayleighWave.C src/MaterialPfile.C
    src/Filter.C src/Polynomial.C src/SecondOrderSection.C src/time_functions.C src/Qspline.C
    src/MaterialIfile.C src/GeographicProjection.C src/Image3D.C src/ESSI3D.C src/ESSI3DHDF5.C
    src/MaterialVolimagefile.C src/MaterialRfile.C src/MaterialSfile.C src/AnisotropicMaterialBlock.C
    src/sacutils.C src/DataPatches.C src/addmemvarforcing2.C src/consintp.C src/oddIoddJinterp.C
    src/evenIoddJinterp.C src/oddIevenJinterp.C src/evenIevenJinterp.C src/CheckPoint.C src/geodyn.C
    src/AllDims.C src/Patch.C src/RandomizedMaterial.C src/MaterialInvtest.C src/sw4-prof.C
    src/sachdf5.C src/readhdf5.C src/TestTwilight.C src/TestPointSource.C src/curvilinear4sgwind.C
    src/TestEcons.C src/GridGenerator.C src/GridGeneratorGeneral.C src/GridGeneratorGaussianHill.C
    src/CurvilinearInterface2.C src/SfileOutput.C src/pseudohess.C src/fastmarching.C src/solveTT.C
    src/rhs4th3point.C src/MaterialGMG.C src/addsgdc.C src/bcfortc.C src/bcfortanisgc.C 
    src/bcfreesurfcurvanic.C src/boundaryOpc.C src/energy4c.C src/checkanisomtrlc.C src/computedtanisoc.C 
    src/curvilinear4sgc.C src/gradientsc.C src/randomfield3dc.C src/innerloop-ani-sgstr-vcc.C 
    src/ilanisocurvc.C src/rhs4curvilinearc.C src/rhs4curvilinearsgc.C src/rhs4th3fortc.C src/solerr3c.C 
    src/testsrcc.C src/rhs4th3windc.C src/tw_aniso_forcec.C src/tw_aniso_force_ttc.C src/velsumc.C 
    src/twilightfortc.C src/twilightsgfortc.C src/tw_ani_stiffc.C src/anisomtrltocurvilinearc.C 
    src/scalar_prodc.C src/updatememvarc.C src/addsg4windc.C src/bndryOpNoGhostc.C src/rhs4th3windc2.C
   )


set(SW4_NATIVE_FORTRAN_SOURCES src/rayleighfort.f src/lamb_exact_numquad.f)

set(SW4_RAJA_SOURCES
  src/EW.C
  src/Sarray.C
  src/version.C
  src/parseInputFile.C
  src/ForcingTwilight.C
  src/curvilinearGrid.C
  src/boundaryOp.f
  src/bndryOpNoGhost.f90
  src/bcfort.f
  src/twilightfort.f
  src/rhs4th3fort.f
  src/parallelStuff.C
  src/Source.C
  src/MaterialProperty.C
  src/MaterialData.C
  src/material.C
  src/setupRun.C
  src/solve.C
  src/solerr3.f
  src/Parallel_IO.C
  src/Image.C
  src/GridPointSource.C
  src/MaterialBlock.C
  src/testsrc.f
  src/TimeSeries.C
  src/sacsubc.C
  src/SuperGrid.C
  src/addsgd.f
  src/velsum.f
  src/rayleighfort.f
  src/energy4.f
  src/TestRayleighWave.C
  src/MaterialPfile.C
  src/Filter.C
  src/Polynomial.C
  src/SecondOrderSection.C
  src/time_functions.C
  src/Qspline.C
  src/lamb_exact_numquad.f
  src/twilightsgfort.f
  src/MaterialIfile.C
  src/MaterialGMG.C
  src/GeographicProjection.C
  src/rhs4curvilinear.f
  src/curvilinear4.f
  src/rhs4curvilinearsg.f
  src/curvilinear4sg.f
  src/gradients.f
  src/Image3D.C
  src/MaterialVolimagefile.C
  src/MaterialRfile.C
  src/randomfield3d.f
  src/innerloop-ani-sgstr-vc.f
  src/bcfortanisg.f
  src/AnisotropicMaterialBlock.C
  src/checkanisomtrl.f
  src/computedtaniso.f
  src/sacutils.C
  src/ilanisocurv.f
  src/anisomtrltocurvilinear.f
  src/bcfreesurfcurvani.f
  src/tw_ani_stiff.f90
  src/tw_aniso_force.f
  src/tw_aniso_force_tt.f
  src/updatememvar.f90
  src/addmemvarforcing2.C
  src/addsg4wind.f90
  src/consintp.C
  src/scalar_prod.f90
  src/oddIoddJinterp.C
  src/evenIoddJinterp.C
  src/oddIevenJinterp.C
  src/evenIevenJinterp.C
  src/CheckPoint.C
  src/Mspace.C
  src/RandomizedMaterial.C
  src/AllDims.C
  src/Patch.C
  src/ESSI3D.C
  src/MaterialSfile.C
  src/MaterialInvtest.C
  src/geodyn.C
  src/ESSI3DHDF5.C
  src/sachdf5.C
  src/readhdf5.C
  src/CurvilinearInterface2.C
  src/TestEcons.C
  src/TestTwilight.C
  src/TestPointSource.C
  src/curvilinear4sgwind.C
  src/GridGeneratorGeneral.C
  src/GridGeneratorGaussianHill.C
  src/GridGenerator.C
  src/RHS43DEV.C
  src/curvilinear4sgcX1.C
  src/curvilinear4sgcSF.C
  src/SfileOutput.C
  src/addsgdc.C
  src/bcfortc.C
  src/bcfortanisgc.C
  src/bcfreesurfcurvanic.C
  src/boundaryOpc.C
  src/energy4c.C
  src/checkanisomtrlc.C
  src/computedtanisoc.C
  src/curvilinear4sgc.C
  src/gradientsc.C
  src/randomfield3dc.C
  src/innerloop-ani-sgstr-vcc.C
  src/ilanisocurvc.C
  src/rhs4curvilinearc.C
  src/rhs4curvilinearsgc.C
  src/rhs4th3fortc.C
  src/solerr3c.C
  src/testsrcc.C
  src/tw_aniso_forcec.C
  src/tw_aniso_force_ttc.C
  src/velsumc.C
  src/twilightfortc.C
  src/twilightsgfortc.C
  src/tw_ani_stiffc.C
  src/anisomtrltocurvilinearc.C
  src/scalar_prodc.C
  src/updatememvarc.C
  src/addsg4windc.C
  src/bndryOpNoGhostc.C
  src/rhs4th3windc2.C
  src/rhs4th3windc.C
)
