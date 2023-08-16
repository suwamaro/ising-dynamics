import numpy as np
import toml
import glob

pfilen = 'parameters.toml'

def calculate_array_average(ifilen, postprocess=None, colidx=0):
    # Extracting parameters
    params = toml.load(open(pfilen))
    L = params['L']
    N = L * L

    # Result
    result = None
    
    # Assume the directory names are 'observables' + numbers
    obsdirs = glob.glob('observables/observables[0-9]*')
    for obsdir in obsdirs:
        df = np.genfromtxt(obsdir+'/'+ifilen, usecols=colidx)
        df = df.reshape((len(df),1))
        if result is None:
            result = df
        else:
            result = np.hstack((result, df))

    # Postprocess
    if postprocess == 'absolute':
        result = np.absolute(result)

    # Taking the average
    mu = np.average(result, axis=1)
    sigma = np.std(result, axis=1, ddof=1) / np.sqrt(len(obsdirs))

    # Time
    nt = len(mu)
    if colidx == 0:
        ts = np.arange(nt)
    else:
        ts = np.genfromtxt(obsdirs[0]+'/'+ifilen, usecols=0)

    # Combining results
    result = np.hstack((ts.reshape((nt,1)), mu.reshape((nt,1)), sigma.reshape((nt,1))))

    # Output
    np.savetxt('observables/'+ifilen, result, fmt='%.10g')



