target_sources(${QUICC_CURRENT_MODEL_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IRTCBackend.cpp
  IRTCModel.cpp
  IRTCState.cpp
  IRTCVisualization.cpp
  Momentum.cpp
  MomentumKernel.cpp
  Transport.cpp
  TransportKernel.cpp
  )

add_subdirectory(Explicit)
add_subdirectory(Implicit)
