include(ConfigureModel)
quicc_add_model(Boussinesq/Sphere/RTC
  TYPES Explicit Implicit
  BROKENMPI Tubular Single1D)
