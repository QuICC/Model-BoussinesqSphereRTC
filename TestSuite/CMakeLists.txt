include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Config"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "State"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Visualization"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
  FILTER "algorithm" "split_equation"
  DATAFILTER "split_equation"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:tubular" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
  FILTER "algorithm" "split_equation"
  DATAFILTER "split_equation"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Config"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:serial"
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "State"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:serial"
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Visualization"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:serial"
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:serial"
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  )

if(TARGET QuICC::Stability)
  quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereRTC"
    TYPE "Stability"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Stability"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
    VARIANTS "framework/parallel/algorithm:tubular"
    MPIRANKS 1
  )
endif()
