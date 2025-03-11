set(type "Explicit")
string(TOLOWER "quicc_${QUICC_CURRENT_MODEL_DIRNAME}" _model_lib)
string(TOLOWER "${_model_lib}_${type}" modLib)
set(ModelId "BoussinesqSphereRTC/${type}")

include(ModelFunctions)
# Config executable
quicc_setup_exe(${ModelId}
  POSTFIX "Config"
  EXESOURCE "WriteConfig.cpp"
  MODELLIB "${modLib}")
