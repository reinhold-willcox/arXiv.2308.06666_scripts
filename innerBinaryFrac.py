import numpy as np
import scipy as sp
import scipy.stats as ss

def getTruncatedPoisson(mu, N=int(1e5)):

    # mu must be np array
    if not isinstance(mu, np.ndarray):
        if isinstance(mu, (list, tuple)):
            mu = np.array(mu)
        else: # scalar, hopefully
            mu = np.array([mu])    
    # mu cannot be non-positive!
    assert np.all(mu>0)

    xx = ss.poisson.rvs(mu=mu, size=(N, mu.shape[0])).T
    trunccts = []
    for xi in xx:
        _, poissonCts = np.unique(xi, return_counts=True)
        truncatedCts = poissonCts[:4]
        trunccts.append(truncatedCts/np.sum(truncatedCts)) # f0, f1, f2, f3
    FN = np.array(trunccts)
    return FN

    

def getPredictedMultiplicityFrequency(mu):
    FN = getTruncatedPoisson(mu)
    all_sums = np.sum(np.array([ii*FN[:,ii] for ii in range(FN.shape[1])]), axis=0)
    return all_sums


def getMuFuncOfMultiplicityFrequency():
    mu = np.linspace(0, 10, 100)[1:]
    fmult = getPredictedMultiplicityFrequency(mu)
    # Artificially add on the max value - these values were chosen to try to continue the fmult(mu) trend naturally
    mu = np.append(mu, 20)
    fmult = np.append(fmult, 3) 
    return sp.interpolate.interp1d(x=fmult, y=mu, fill_value='extrapolate', kind='slinear') 


def getMultFreqFuncOfM1FromObs(m1_obs, fmult_obs):
    return sp.interpolate.interp1d(x=m1_obs, y=fmult_obs, fill_value='extrapolate', kind='quadratic')


def buildTrueBinaryFuncOfM1(fN):
    trueBinFracVals = fN[1]/(np.sum(fN[1:], axis=0)) # True bin frac is the number of true bins out of all binaries + triples + higher
    return sp.interpolate.interp1d(M1, trueBinFracVals, fill_value='extrapolate')


def trueBinaryFuncOfM1(M1, Mu):
    # Get probability of N companions
    fN = getTruncatedPoisson(Mu).T
    # Sum over to get true binary fraction per M1
    trueBinFracVals = fN[1]/(np.sum(fN[1:], axis=0)) # True bin frac is the number of true bins out of all binaries + triples + higher
    return sp.interpolate.interp1d(M1, trueBinFracVals, fill_value='extrapolate')


def getTrueBinaryFractionFuncOfM1():

    fmult_obs = [.5, .84, 1.3, 1.6,   2.0, 2.4]
    m1_obs =    [1,  3.5, 7,   12.5,  30, 50]    
    
    # Calculate the functions 
    multFreqFuncOfM1 = getMultFreqFuncOfM1FromObs(m1_obs, fmult_obs)
    muFuncOfMultFreq = getMuFuncOfMultiplicityFrequency()

    # Get a dummy array for the primary mass
    M1=np.logspace(np.log10(.5), np.log10(100), 100)
    # Combine the functions to get
    Fm = multFreqFuncOfM1(M1)
    Mu = muFuncOfMultFreq(Fm)
            
    return trueBinaryFuncOfM1(M1, Mu)





