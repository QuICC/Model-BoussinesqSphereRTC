list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

# Set Model directory name
get_filename_component(QUICC_CURRENT_MODEL_DIRNAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(QUICC_CURRENT_MODEL_DIR "Models/${QUICC_CURRENT_MODEL_DIRNAME}")
string(TOLOWER "quicc_${QUICC_CURRENT_MODEL_DIRNAME}" QUICC_CURRENT_MODEL_LIB)

include(ConfigureModel)

quicc_add_model(Boussinesq/Sphere/RTC
  TYPES Explicit Implicit
  MODEL_DIRNAME ${QUICC_CURRENT_MODEL_DIRNAME})

if(TARGET QuICC::Tools::Simulation)
  list(APPEND CMAKE_MODULE_PATH "${QUICC_TOOL_CMAKE_MODULE_PATH_Simulation}")
  message(DEBUG "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  include(ConfigureSimulation)
  quicc_add_simulation(Boussinesq/Sphere/RTC
    TYPES Explicit Implicit
    MODEL_DIRNAME ${QUICC_CURRENT_MODEL_DIRNAME})
endif()

if(TARGET QuICC::Tools::Stability)
  list(APPEND CMAKE_MODULE_PATH "${QUICC_TOOL_CMAKE_MODULE_PATH_Stability}")
  message(DEBUG "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  include(ConfigureStability)
  quicc_add_stability(Boussinesq/Sphere/RTC
    TYPES Implicit
    MODEL_DIRNAME ${QUICC_CURRENT_MODEL_DIRNAME})
endif()

# Teststuite, not all models have a test
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d/setup/TestSuite.cmake")
  include(setup/TestSuite)
endif()
