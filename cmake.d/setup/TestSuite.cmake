option(QUICC_TESTSUITE_BOUSSINESQSPHERERTC "Enable BoussinesqSphereRTC model testsuite?" OFF)
if(QUICC_TESTSUITE_BOUSSINESQSPHERERTC)
  add_subdirectory(TestSuite)
endif()
