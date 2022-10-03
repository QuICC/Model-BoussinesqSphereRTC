include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  RUNID "_single1d"
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  RUNID "_single2d"
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  RUNID "_tubular"
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  RUNID "_single1d"
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  RUNID "_single2d"
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  RUNID "_tubular"
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  )
