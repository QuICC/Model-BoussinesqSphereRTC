set(QUICC_WORK_DIR
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks")

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/Models/_refdata")

set(_impl "Explicit")
set(_benchmarkExe "BoussinesqSphereRTC${_impl}Model")
set(_benchmarkTarget "Benchmark${_benchmarkExe}")

set(_refdir "${QUICC_WORK_DIR}/_refdata/${_impl}")
set(_rundir "${QUICC_WORK_DIR}/_data/${_impl}")
add_custom_target(${_benchmarkTarget} ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink
    "${_refdir}/parameters.cfg"
    "${_rundir}/parameters.cfg"
  COMMAND ${CMAKE_COMMAND} -E create_symlink
    "${_refdir}/state_initial.hdf5"
    "${_rundir}/state_initial.hdf5"
  )
add_dependencies(${_benchmarkTarget} ${_benchmarkExe})

# Fetch reference data
include(FetchBenchmarkReference)
quicc_fetch_benchmark_reference(
  ${_benchmarkTarget}
  MODEL BoussinesqSphereRTC
  FILENAME "${_impl}.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

set(_testname "Benchmark${_benchmarkExe}")
add_test(
  NAME ${_benchmarkTarget}
  COMMAND "${_benchmarkExe}"
  WORKING_DIRECTORY
  "${_rundir}"
  )
