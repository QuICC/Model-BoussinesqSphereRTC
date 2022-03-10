import numpy as np

ref_dir = "../Explicit/"
data_dir = "./"
results = np.zeros(2, dtype='i8')

eps = 1e-12

def tableTest(fname, usecols = None):
    # Validate nusselt number
    print(f'Validating {fname}')
    ref = np.genfromtxt(ref_dir + fname, usecols=usecols)
    data = np.genfromtxt(data_dir + fname, usecols=usecols)
    res = np.zeros(2, dtype='i8')

    if(ref.shape == data.shape):
        status = "passed"
    else:
        status = "failed"
        res[1] += 1
    res[0] += 1
    print(f'\tChecking size: {status}')

    max_err = np.max(np.abs(ref - data))

    if(max_err < eps):
        status = "passed"
    else:
        status = "failed"
        res[1] += 1
    res[0] += 1
    print(f'\tMaximum error ({max_err}): {status}')

    return res

# Temperature
#   Energy
results += tableTest("temperature_energy.dat")
#   L spectrum
results += tableTest("temperature_l_spectrum0000.dat")
results += tableTest("temperature_l_spectrum0100.dat")
#   M spectrum
results += tableTest("temperature_m_spectrum0000.dat")
results += tableTest("temperature_m_spectrum0100.dat")
#   M spectrum
results += tableTest("temperature_n_spectrum0000.dat")
results += tableTest("temperature_n_spectrum0100.dat")
# Kinetic
#   energy
results += tableTest("kinetic_energy.dat")
#   L spectrum
results += tableTest("kinetic_l_spectrum0000.dat")
results += tableTest("kinetic_l_spectrum0100.dat")
#   M spectrum
results += tableTest("kinetic_m_spectrum0000.dat")
results += tableTest("kinetic_m_spectrum0100.dat")
#   N spectrum
results += tableTest("kinetic_n_spectrum0000.dat")
results += tableTest("kinetic_n_spectrum0100.dat")
# Nusselt number
results += tableTest("nusselt.dat")
# Angular momentum
results += tableTest("angular_momentum.dat")
# CFL
results += tableTest("cfl.dat", usecols=(0,1,3,5,6,7,8,9))

print("")
if(results[1] == 0):
    print(f'All tests passed')
else:
    t = 'test'
    if results[1] > 1:
        t += 's'
    print(f'{results[1]} {t} failed out of {results[0]}')
