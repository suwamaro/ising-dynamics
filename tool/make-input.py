import sys
import os
import shutil
import os.path
import numpy as np
from copy import deepcopy

CURRENT = os.getcwd()
HOME = CURRENT
print('HOME =', HOME)

# Parameters
DATE = '-1202'
SUFF = ''
N0 = 100
nsim = 10
nskip = 1
Ns = np.arange(N0, N0+nsim*nskip, step=nskip)
var_name = 'kT'
dp = 0.1
parami = 1.0
paramf = parami + dp * (nsim-0.5)
vars = np.arange(parami, paramf, dp)
nd = 8  # precision for the parameters
nf = '{:.'+str(nd)+'g}'

params = {
    'L': 32,
    'J1': 1.0,
    'J2': 1.5,
    'J3': 0.75,
    'Jp': -1.0,
    'kT': 2.04,
    'n_mcs': 1024,
    'fast_code': 'true',
    'measure_energy': 'false',
    'dump_config': 'true',
    'n_dump': 8,
    'seed': 0,
    'n_replicas': 64,
}

print('Ns:', Ns)
print('var_name:', var_name)
print('vars:', vars)
if len(Ns) != len(vars):
    print('Ns and vars have different lengths.')
    print('len(Ns) =', len(Ns), '; len(vars) =', len(vars))
    quit()
print('Number of parameter sets:', len(vars))

ofilen = 'parameters.toml'
OVERWRITE = False

argv = sys.argv
argv.pop(0)

# Usage
usage = 'Usage: python {} directory numbers.'\
        .format(__file__)

# Options
options = [option for option in argv if option.startswith('-')]
if '--help' in options:
    exit(usage)

# Checking if directories exist.
Nss = deepcopy(Ns)
Ns = []
for n in Nss:
    simn = HOME+'/sim' + DATE + '-' + str(n) + SUFF
    if os.path.exists(simn):
        Ns.append(n)

if not OVERWRITE and Ns:
    for n in Ns:
        exit(HOME+'/sim'+str(DATE)+'-'+str(n)+SUFF+' already exists.')
else:
    Ns = deepcopy(Nss)

def write_params():
    out = open(ofilen, 'w')
    for key, value in params.items():
        out.write(key + ' = ' + str(value) + '\n')
    out.close()

for i,n in enumerate(Ns):
    simdir =  HOME + '/sim' + str(DATE) + '-' + str(n) + SUFF
    print(simdir)
    os.mkdir(simdir)
    os.chdir(simdir)
    params['seed'] = str(int(DATE[1:]))+str(n)
    params[var_name] = str(nf.format(vars[i]))
    write_params()
    os.chdir('..')
