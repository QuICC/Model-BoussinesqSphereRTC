set(_vel_basis "VALUE_TOR_VALUE_POL" "VALUE_TOR_INSULATING_POL" "VALUE_TOR_NS_POL")

quicc_create_option(NAME QUICC_BESSEL_VELOCITY_BC
                    OPTS ${_vel_basis}
                    LABEL "Bessel basis for velocity"
                    ADVANCED)
quicc_add_definition(QUICC_BESSEL_VELOCITY_BC)
