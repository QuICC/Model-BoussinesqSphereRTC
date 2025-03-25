include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
  FILTER "algorithm" "split_equation"
  DATAFILTER "split_equation"
  )

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
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:tubular" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
  FILTER "algorithm" "split_equation"
  DATAFILTER "split_equation"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:serial"
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
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  )

if(TARGET QuICC::Stability)
  include(ConfigureStabilityBenchmark)
  quicc_add_stability_benchmark(Implicit
    MODEL "BoussinesqSphereRTC"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Stability"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
    VARIANTS "framework/parallel/algorithm:tubular"
    MPIRANKS 1
    MODEL_DIRNAME ${QUICC_CURRENT_MODEL_DIRNAME}
  )
endif()
