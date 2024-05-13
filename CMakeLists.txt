include (cmake.d/setup/BesselBasis.cmake)

include(ConfigureModel)
quicc_add_model(Boussinesq/Sphere/RTC
  TYPES Explicit Implicit)
