# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Want to get distributions of the Moe & Distefano parameters as functions of singles binaries and inner binaries of triples
#
#

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.stats as ss

# +
 
### Constants
AU_per_m = (6.68459e-12)
s_per_dy = (3600*24)
dy_per_yr = (365.2422)
s_per_yr = s_per_dy*dy_per_yr
kg_per_Msun = (1.989e30)

G = 6.67430e-11 # [m^3 /s^2 /kg]
G_au_dy_msun = G * AU_per_m**3 * s_per_dy**2 * kg_per_Msun # [AU^3 /dy^2 /Msun]

#from astroUtils import PtoA
def PtoA(P, Mtot): 
     """ 
     method to convert period to semi-major-axis 
     IN:
         P : period in days
         Mtot: total mass in Msun 
     OUT:
         semi-major axis in AU
     """
     const = G_au_dy_msun/(4*(np.pi)**2)
     return np.cbrt(const*Mtot*np.square(P)) # AU



# -

# ## Estimate the fraction of binaries from Moe & DiStefano

# +

def getTrucatedPoisson(mu, N=int(1e5)):
    
    # Generate Poisson RVs and get the relative frequencies of bins [0,3]
    xx = ss.poisson.rvs(mu=mu, size=N)
    _, poissonCts = np.unique(xx, return_counts=True)
    truncatedCts = poissonCts[:4]
    return truncatedCts/np.sum(truncatedCts) # f0, f1, f2, f3

def getPredictedMultiplicityFrequency(mu):
    try:
        ff = getTrucatedPoisson(mu)
        return sum([ii*ff[ii] for ii in range(len(ff))])
    except:
        all_sums = []
        for m in mu:
            ff = getTrucatedPoisson(m)
            all_sums.append(sum([ii*ff[ii] for ii in range(len(ff))]))
        return all_sums
    
def getMultiplicityFractionArray():
    
    mu = np.linspace(.1, 10, 100)
    fm = getPredictedMultiplicityFrequency(mu)
    mu_func = sp.interpolate.interp1d(x=fm, y=mu, fill_value=(0, 0)) 

    # want to interpolate to get a function that gives mu as a function of fmult
    fmult_vals = [.5, .84, 1.3, 1.6,   2.0, 2.4]
    m1_vals =    [1,  3.5, 7,   12.5,  30, 50]

    mu_vals = mu_func(fmult_vals)
    fm_vals = np.linspace(.1, 2.5, 100)

    fmult_func_of_m1 = sp.interpolate.interp1d(x=m1_vals, y=fmult_vals, fill_value='extrapolate')
    fmult_func_of_m1_2 = sp.interpolate.interp1d(x=m1_vals, y=fmult_vals, fill_value='extrapolate', kind='quadratic')

    # Combine to get Fn dist as a function of m1
    M1=np.logspace(np.log10(.5), np.log10(50), 100)
    fmult_vals = fmult_func_of_m1_2(M1)
    mu_vals = mu_func(fmult_vals)

    fN = []
    for m1, mu in zip(M1, mu_vals):
        try:
            fn_vals = getTrucatedPoisson(mu)
        except:
            continue        
        fN.append(fn_vals)

    fN = np.array(fN).T
    
    return fN, M1


def getMultiplicityFunctionsOfMass():
    fN, M1 = getMultiplicityFractionArray()
    fN_func = []
    for ii, fn in enumerate(fN):
        fN_func.append(sp.interpolate.interp1d(x=M1, y=fn, fill_value='extrapolate'))
    return M1, fN_func


def buildTrueBinaryFunction():
    fN, M1 = getMultiplicityFractionArray()
    trueBinFrac = fN[1]/(np.sum(fN[1:], axis=0))
    calculateTrueBinFracAsFunctionOfMass = sp.interpolate.interp1d(M1, trueBinFrac, bounds_error=False, fill_value=(.8, .1))
    
    return calculateTrueBinFracAsFunctionOfMass


# -
calculateTrueBinFracAsFunctionOfMass = buildTrueBinaryFunction()


# +
# Plot above results, ensure they match

fig, axes = plt.subplots(ncols=2, figsize=(20, 8))
fs = 24
colors = ['red', 'g', 'b', 'purple']

ax = axes[0]
M1, fN_func = getMultiplicityFunctionsOfMass()
for ii, fn in enumerate(fN_func):
    ax.plot(M1, fn(M1), color=colors[ii])
ax.set_xscale('log')
ax.set_xlabel(r'Primary Mass $M_1 (M_\odot)$', fontsize=fs)
ax.set_ylabel(r'$N_\mathrm{comp}$', fontsize=fs)
ax.set_title('Number of companions per primary mass')



ax = axes[1]
m1 = np.logspace(0, 2, 101)
fbin = calculateTrueBinFracAsFunctionOfMass(m1)
ax.plot(m1, fbin)
ax.set_xscale('log')
ax.set_ylim(0, 1)
ax.set_xlabel(r'Primary Mass $M_1 (M_\odot)$', fontsize=fs)
ax.set_ylabel(r'$P_{{bin}}(M_1)$', fontsize=fs)
ax.set_title('Probability that a sampled binary is a true binary')


# -








# +
with open('grid_moedistefano_nSamples100000.txt', 'r') as f:
    lines = f.readlines()
    lines = [line.rstrip() for line in lines]

M1 = []
M2 = []
P = []
E = []

for line in lines:
    _, m1, _, m2, _, p, _, e = line.split()
    M1.append(m1)
    M2.append(m2)
    P.append(p)
    E.append(e)

def convertToArray(datalist):
    return np.array(datalist, dtype=float)

M1 = convertToArray(M1)
M2 = convertToArray(M2)
P = convertToArray(P) # days
E = convertToArray(E)

Q = M2/M1
A = PtoA(Mtot=M1+M2, P=P) # AU
RP = A*(1-E)

# Get fractions of singles and binaries
maskSingles = P == 100000
maskBinaries = ~maskSingles

# Get probability of inner binary
probTrueBinary = calculateTrueBinFracAsFunctionOfMass(M1)
probTrueBinary[maskSingles] = 0 # singles are not true binaries

randDraw = np.random.rand(len(M1))
maskTrueBinaries = probTrueBinary > randDraw
maskInnerBinaries = maskBinaries & ~maskTrueBinaries



# +
# Get histograms

fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(20, 20))
axes = axes.flatten()

masks = [maskSingles, maskTrueBinaries, maskInnerBinaries]
colors = ['deepskyblue', 'forestgreen', 'orange']
fs_label=50
fs_nums = 30

for ax in axes:
    ax.tick_params(axis='both', which='major', labelsize=fs_nums)


# M1
ax = axes[0]
ax.hist([np.log10(M1)[mask] for mask in masks], histtype='barstacked', density=True, bins=np.linspace(0.5, 2, 41), color=colors)
ax.set_xlabel(r'Primary Mass $\log(M_1/M_\odot)$', fontsize=fs_label)
ax.set_ylabel('Sample density', fontsize=fs_label)
ax.set_xlim(.6, 2.05)

# q
ax = axes[1]
ax.hist([Q[mask] for mask in masks[1:]], histtype='barstacked', density=True, bins=np.linspace(0, 1, 21), color=colors[1:])
ax.set_xlabel(r'Mass Ratio $q = M_2/M_1$', fontsize=fs_label)
#ax.set_ylabel('Sample density', fontsize=fs_label)

# rP
ax = axes[2]
ax.hist([np.log10(RP)[mask] for mask in masks[1:]], histtype='barstacked', density=True, bins=np.linspace(-2, 4, 41), color=colors[1:])
ax.set_xlabel(r'Periapsis $\log(r_P/\mathrm{AU})$', fontsize=fs_label)
ax.set_ylabel('Sample density', fontsize=fs_label)
ax.set_xlim(-1.8, 4.5)

# pie
ax = axes[3]
"""
labels = [r'Singles', r'True Binaries', 'Inner Binaries', 'Hierarchical Systems']
textprops={'fontsize':fs_label, 'color':'k', 'weight':'bold'}
fracs = [np.sum(mask) for mask in masks]
ax.pie(fracs, colors=colors,  radius=1.3)#, labels=labels, labeldistance=.3, )
ax.text(s=labels[0], x=.4, y=.3, **textprops)
ax.text(s=labels[1], x=-1.2, y=.2, **textprops)
ax.text(s=labels[2], x=-.7, y=-.4, **textprops)
#ax.text(s=labels[3], x=-.9, y=-.6, **textprops)
"""
fracs = [np.sum(mask) for mask in masks]
fracs /= np.sum(fracs)
labels = [r'Singles', r'True Binaries', 'Inner Binaries', 'Hierarchical Systems']
fstr = ['{:.2f}\%'.format(frac*100) for frac in fracs]

textprops={'fontsize':fs_label, 'color':'k', 'weight':'bold'}
ax.pie(fracs, colors=colors,  radius=1.3)#, labels=labels, labeldistance=.3, )
ax.text(s=labels[0], x=.4, y=.1, **textprops)
ax.text(s=fstr[0], x=.4, y=.3, **textprops)
ax.text(s=labels[1], x=-1.2, y=.2, **textprops)
ax.text(s=fstr[1], x=-.8, y=.4, **textprops)
ax.text(s=labels[2], x=-.7, y=-.6, **textprops)
ax.text(s=fstr[2], x=-.3, y=-.4, **textprops)

#ax.text(s=labels[3], x=-.9, y=-.6, **textprops)


# full title
#fig.suptitle('Initial Distributions from Moe \& Di Stefano (2017)', fontsize=30)

fig.tight_layout(pad=5)

#plt.savefig('../plots/final_plots/moedistef_ICs.pdf', format='pdf', bbox_inches='tight')
plt.show()


# -
#

