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
  VARIANTS "framework/parallel/algorithm:serial"
  FILTER "algorithm"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "setup/boundary/scheme:Galerkin"
  FILTER "algorithm" "scheme"
  VALIDATIONFILTER "boundary/scheme"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2"
  FILTER "algorithm" "scheme"
  DATAFILTER "scheme"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/boundary/scheme:Galerkin"
  FILTER "algorithm" "scheme"
  DATAFILTER "timestepping/scheme"
  VALIDATIONFILTER "boundary/scheme"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
  FILTER "algorithm" "timestepping/scheme" "split_equation"
  DATAFILTER "timestepping/scheme"
  VALIDATIONFILTER "split_equation"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Model"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_RK_CB3b"
  FILTER "algorithm" "scheme"
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
  FILTER "algorithm" "timestepping/scheme" "split_equation"
  DATAFILTER "timestepping/scheme"
  VALIDATIONFILTER "split_equation"
  MPIRANKS 6
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Config"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "State"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqSphereRTC"
  TYPE "Visualization"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
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
    TYPE "StabilityConfig"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Stability"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:tubular"
    MPIRANKS 1
  )

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
