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

import astropy.io.fits as fits

fnames = ['../data/ge2015_tab3/J_ApJ_812_40_table3.dat.fits.gz', '../data/ge2020_tab3/J_ApJ_899_132_table3.dat.fits.gz']


def getAllData(fname):
    with fits.open(fname) as hdul:
        hdul.info()
        print()
        allData = hdul[1].data
        print(repr(hdul[1].header))
    return allData



allDatas = [getAllData(fname) for fname in fnames]

# visualize the data
# for ii in range(20):
#    print(allData[ii])


# Print a COMPAS style table of the mass radius qadic relations
def getTable(allData):
    try:
        logRad = allData['logRad']
    except:
        logRad = allData['Rad']    
    table = np.array([allData['Mass'], logRad, allData['qad'], allData['qadic'], allData['zetaad'], allData['zetaadic']])
    print(table.shape)
    return table



tables = [getTable(allData) for allData in allDatas]





uniqueMasses = np.unique(tables[0][0])
largeMasses = uniqueMasses[uniqueMasses>3.9]


def getDataDict(table):
    dataDict = dict()

    dependent_params = [ 'logRad', 'qad', 'qadic', 'zetaad', 'zetaadic' ]

    for mm in uniqueMasses:
        mask = table[0] == mm    
        val = np.array([table[jj+1][mask] for jj in range(len(dependent_params))])
        dataDict.update({mm: val})
    
    return dataDict


dataDicts = [getDataDict(table) for table in tables]


# +

def zeta_SPH(m):
    # m is m_core/m_total
    
    zeta_sph1 = (2/3)*(m/(1-m)) - (1/3)*((1-m)/(1+2*m))
    
    one_minus_m_neg6 = np.power(1-m, -6)    
    rhs = -0.03*m + 0.2* (m/(1+one_minus_m_neg6))

    return zeta_sph1 + rhs

def qClaeysGB(m):
    # m is mcore/m_total
    
    m5 = m*m*m*m*m
    return 2.13/(1.67-.3+2*m5)

# +
# For specific masses, plot the various params as a function of log radius


def plotZetaQad(ax, mm, dataDict, plotQs=True):
    logR, qad, qadic, zetaad, zetaadic = dataDict[mm]

    #yvals = [zetaad, qad, zetaadic, qadic]
    yvals = [zetaad, 1.0/qad, zetaadic, 1.0/qadic]
    linestyles = [':', ':', '-.', '-.']
    alpha = [1, 1, 0.5, 0.5]
    col1, col2 = 'blue', 'blue'
    colors = [col1, col1, col2, col2]
    labels = [r'$\zeta$', r'$q$', r'$\tilde{\zeta}$', r'$\tilde{q}$']

    # Remove any R values that are less than previous values
    mask = ~(logR < np.maximum.accumulate(logR))
    
    #for ii in [0, 2, 1, 3]:
    if plotQs:
        myInds = [1, 3]
    else:
        myInds = [0, 2]
        
    for ii in myInds:
        yval = yvals[ii]
        if plotQs:
            mask1 = mask & (yval < 1) & (yval > 0)    
        else:
            #yval[yval>20] = 20
            mask1 = mask
        ax.plot(logR[mask1], yval[mask1], ls=linestyles[ii], color=colors[ii], alpha=alpha[0], label=labels[ii], lw=2)
    
    if False: #plotQs:
        if mm < 12:
            loc = 'upper left'
        else:
            loc = 'upper right'
        ax.legend(fontsize=18, loc=loc, framealpha=1)

    ax.set_title(r'{}$M_\odot$'.format(mm))
    


# +
def convertQCritToZetaAtBeta1(q):
    first = -2*(1-1/q)
    q13 = np.cbrt(q)
    q23 = q13*q13
    sec_num = 1.2 + 1/(1/q13 + 1/q23)
    sec_den = 0.6 + q23*np.log(1+1/q13)
    second = (1/3)*(2-sec_num/sec_den)*(1/q + 1)
    return first + second
    
#testing
#q = np.linspace(0.1, 1, 10)
#xout = convertQCritToZetaAtBeta1(q)
#xout


# +
def get_h5(ii):
    fname = '../data/single_star_grid_output//Detailed_Output/BSE_Detailed_Output_{}.h5'.format(ii)
    return h5.File(fname, 'r')


def plotTransitionRadii(ax, mm, bounds=None, addRadialLines=False, addZetaSPH=False, addQCritClaeys=False, addZetaClaeys=False):
    for ii in range(100):
        myf = get_h5(ii)
        m0 = myf['Mass(1)'][()][0]
        if mm == m0:
            break

    M = myf['Mass(1)'][()]
    Mc = myf['Mass_He_Core(1)'][()]
    R = myf['Radius(1)'][()]
    mask1 = ~(R < np.maximum.accumulate(R))
    stype = myf['Stellar_Type(1)'][()]
    stype_upper = 5 if mm < 12 else 7
    mask2 = np.append(np.array([False]), stype[mask1][1:]>stype[mask1][:-1]) & (stype[mask1]<stype_upper)
    Rads = R[mask1][mask2]
    STs = stype[mask1][mask2]
    logRs = [np.log10(Rad) for Rad in Rads]
    
    
    #if not addZetaSPH:
    if addRadialLines:
        st_str = ['-', '-', 'HG', 'BGB', 'CHeB', 'EAGB', 'TPAGB']
        xcheb = -.15 if mm < 12 else -.3 if mm != 25 else -.55
        x_off = ['-', '-', -.1, -.3, xcheb, -.4, -.1]
        if not addZetaSPH:
            yupp, ylow, yoff = 3/24, 23/24, -.9
            #ypt = -.075

        else:
            yupp, ylow, yoff = 3/25, 1, -0.9
            #ypt = -.075

        for st, logR in zip(STs, logRs):
            ax.axvline(logR, yupp, ylow, color='k', linestyle='--', alpha=0.5)
            ax.text(logR+x_off[st], yoff, s=st_str[st], fontsize=18)#, transform=ax.transAxes)
        
    if addZetaSPH:
        logR = np.log10(R)
        mask = ~(logR < np.maximum.accumulate(logR)) & (logR>bounds[0]) & (logR<bounds[1])
        logR = logR[mask]
        zetaComp = 2*np.ones_like(logR)
        mask1 = (stype > 2)[mask]
        zetaComp[mask1] = zeta_SPH(Mc[mask][mask1]/M[mask][mask1])
        mask2 = (stype == 2)[mask]
        zetaComp[mask2] = 6.5
        # Make transition off HG vertical
        indInsert = np.where(zetaComp == 6.5)[0][-1] + 1 # insert just after the HG stuff
        logR = np.insert(logR, indInsert, logR[indInsert])
        zetaComp = np.insert(zetaComp, indInsert, 6.5)
        # Make transition onto HG vertical
        indInsert = np.where(zetaComp == 6.5)[0][0] # insert just before the HG stuff
        logR = np.insert(logR, indInsert, logR[indInsert])
        zetaComp = np.insert(zetaComp, indInsert, 2)
        ax.plot(logR, zetaComp, 'k-', label=r'$\zeta_{C}$', lw=2)
        #print(zetaComp)
        
    if (addQCritClaeys or addZetaClaeys):
        logR = np.log10(R)
        mask = ~(logR < np.maximum.accumulate(logR)) & (logR>bounds[0]) & (logR<bounds[1])
        logR = logR[mask]
        qComp = 0.625*np.ones_like(logR)
        mask1 = (stype > 2)[mask]
        qComp[mask1] = qClaeysGB(Mc[mask][mask1]/M[mask][mask1])
        mask2 = (stype == 2)[mask]
        qComp[mask2] = .25
        # Make transition off HG vertical
        indInsert = np.where(qComp == .25)[0][-1] + 1 # insert just after the HG stuff
        logR = np.insert(logR, indInsert, logR[indInsert])
        qComp = np.insert(qComp, indInsert, .25)
        # Make transition onto HG vertical
        indInsert = np.where(qComp == .25)[0][0] # insert just before the HG stuff
        logR = np.insert(logR, indInsert, logR[indInsert])
        qComp = np.insert(qComp, indInsert, .625)
        
        if addQCritClaeys:
            ax.plot(logR, qComp, 'r:', label=r'$q_{C}$', lw=2)
        else:
            # Need to convert from q to zeta
            zetaClaeys = convertQCritToZetaAtBeta1(qComp)
            #print(zetaClaeys)
            ax.plot(logR, zetaClaeys, '--', color='darkorange', label='zClaeys')


# +
#uniqueLargeMasses = np.unique(table[0])
fig, axes = plt.subplots(ncols = 5, nrows=3, figsize=(20, 12))


# Plot ZETAs
qMinMax = []
for ii, mm in enumerate(largeMasses):
    ax = axes.flatten()[ii]
    #plotZetaQad(ax, mm, dataDicts[0], '-')
    plotZetaQad(ax, mm, dataDicts[1], plotQs=False)
    xlims = ax.get_xlim()
    ax.set_xlim(*xlims)
    
    fs = 30
    if ii%5 == 0:
        ax.set_ylabel(r'$\zeta_*$', fontsize=fs)
    if ii > 9:
        ax.set_xlabel(r'log($R/R_\odot$)', fontsize=fs)

    #ax.set_ylim(-5, 20)

    rect = mpl.patches.Rectangle(xy=(xlims[0],-4.99), width=xlims[1]-xlims[0], height=4.9, color='w', alpha=1.0)
    ax.add_patch(rect)
    
    plotTransitionRadii(ax, mm, bounds=ax.get_xlim(), addZetaSPH=True, addZetaClaeys=True)
    
    if mm < 24:
        loc = 'upper left'
    else:
        loc = 'upper right'
    #ax.legend(fontsize=18, loc=loc, framealpha=1)

    ax.grid('True', alpha=0.3)
    ax.set_axisbelow(True)



fig.tight_layout()

#plt.savefig('../plots/final_plots/zetas_vs_R.pdf', format='pdf', bbox_inches='tight')


# +
print('--> Use this one! <--')

#uniqueLargeMasses = np.unique(table[0])
fig, axes = plt.subplots(ncols = 5, nrows=3, figsize=(20, 12))


# Plot ZETAs
qMinMax = []
for ii, mm in enumerate(largeMasses):
    ax = axes.flatten()[ii]
    #plotZetaQad(ax, mm, dataDicts[0], '-')
    plotZetaQad(ax, mm, dataDicts[1], plotQs=False)
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.set_xlim(*xlims)
    ax.set_yscale('symlog', linthresh=1)
    ax.set_yticks([-0.5, 0, 0.5, 1, 10, 100, 1000, 10000])
    ax.set_yticklabels(['-0.5', '0', '0.5', '1', '10', r'$10^{2}$', r'$10^{3}$', r'$10^{4}$'])
    minor_ticks = np.concatenate((
        np.linspace(-1, 1, 21),
        np.linspace(1, 10, 11),
        np.linspace(10, 100, 11),
        np.linspace(100, 1000, 11),
        np.linspace(1000, 10000, 11)
    ))
    ax.set_yticks(minor_ticks, minor=True)
    ax.set_ylim(-1, max(100, ylims[1]*2))
    ylims = ax.get_ylim()
    
    fs = 30
    if ii%5 == 0:
        ax.set_ylabel(r'$\zeta_{ad}$', fontsize=fs)
    if ii > 9:
        ax.set_xlabel(r'log($R/R_\odot$)', fontsize=fs)

    #ax.set_ylim(-5, 20)

    rect = mpl.patches.Rectangle(xy=(xlims[0],-1), width=xlims[1]-xlims[0], height=.47, color='w', alpha=1.0)
    ax.add_patch(rect)
    
    rect = mpl.patches.Rectangle(xy=(xlims[0],1), width=xlims[1]-xlims[0], height=ylims[1]-1, color='k', alpha=0.03)
    ax.add_patch(rect)
    
    plotTransitionRadii(ax, mm, bounds=ax.get_xlim(), addZetaSPH=True, addZetaClaeys=True, addRadialLines=True)
    
    ax.plot(xlims, [1, 1], 'k-', alpha=0.25) 
    
    #if mm < 24:
    #    loc = 'upper left'
    #else:
    #    loc = 'upper right'
    #ax.legend(fontsize=18, loc=loc, framealpha=1)

    ax.grid('True', alpha=0.3)
    ax.set_axisbelow(True)



fig.tight_layout()

#plt.savefig('../plots/final_plots/zetas_vs_R.pdf', format='pdf', bbox_inches='tight')
# -








# +
#uniqueLargeMasses = np.unique(table[0])
fig, axes = plt.subplots(ncols = 5, nrows=3, figsize=(20, 12))

# plot Qs

qMinMax = []
for ii, mm in enumerate(largeMasses):
    ax = axes.flatten()[ii]
    #plotZetaQad(ax, mm, dataDicts[0], '-')
    plotZetaQad(ax, mm, dataDicts[1])
    xlims = ax.get_xlim()
    ax.set_xlim(*xlims)
    
    fs = 30
    if ii%5 == 0:
        ax.set_ylabel(r'$q_c$', fontsize=fs)
    if ii > 9:
        ax.set_xlabel(r'log($R/R_\odot$)', fontsize=fs)
    #ax.set_ylim(-0.15, 1.05)
    ax.set_ylim(-0.2, 1.6)
    

    #ax.plot(logR, ypt+.01, 'wo', ms=20)   
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    rect = mpl.patches.Rectangle(xy=(xlims[0],ylims[0]), width=xlims[1]-xlims[0], height=-ylims[0]-.01, color='w', alpha=1.0)
    ax.add_patch(rect)
    
    plotTransitionRadii(ax, mm, bounds=ax.get_xlim(), addQCritClaeys=True)
    
    ax.grid('True', alpha=0.3)
    ax.set_axisbelow(True)



fig.tight_layout()

#plt.savefig('../plots/final_plots/qCrit_ge20.pdf', format='pdf', bbox_inches='tight')
# -







# +
# print line of m vals.
table = tables[1]

firstline = "{" + ", ".join(r"{}".format(m) for m in np.unique(table[0])) + "},"
print(firstline)
print("{")

#print()

for m in np.unique(table[0]):

    mask = table[0] == m
    
    logR = table[1][mask]
    # Only want values on the HG, after the dip at the end of the MS
    mask1 = ~(logR < np.maximum.accumulate(logR))
    
    rr = table[1][mask][mask1]
    qq1 = table[2][mask][mask1]
    qq2 = table[3][mask][mask1]
    zz1 = table[4][mask][mask1]
    zz2 = table[5][mask][mask1]
    
    
    line = "  {{" + ", ".join(r"{}".format(r) for r     in rr) + "}"
    line+= ",  {" + ", ".join(r"{:.3f}".format(1.0/q) for q in qq1) + "}"
    line+= ",  {" + ", ".join(r"{:.3f}".format(1.0/q) for q in qq2) + "}"
    line+= ",  {" + ", ".join(r"{}".format(z) for z     in zz1) + "}"
    line+= ",  {" + ", ".join(r"{}".format(z) for z     in zz2) + "}},"

    fullline = line.replace('<', '{').replace('>', '}')
    print(fullline)
    
    

print("}")
# -










# +
# #ls ../data/single_star_grid_output/Detailed_Output/BS
# -

ii = 90
fname = '../data/single_star_grid_output//Detailed_Output/BSE_Detailed_Output_{}.h5'.format(ii)
myf = h5.File(fname, 'r')
myf.keys()


# +
M = myf['Mass(1)'][()] # 12 Msol star
print(M[0])
Mc = myf['Mass_He_Core(1)'][()]
t = myf['Time'][()]
R = myf['Radius(1)'][()]
stype = myf['Stellar_Type(1)'][()]

fig, axes = plt.subplots(ncols=4, figsize=(25, 8))

param = [M, Mc, R, stype]
name = ['M', 'Mc', 'R', 'St']

for ii in range(4):
    
    ax = axes[ii]
    ax.plot(t, param[ii], label=name[ii])
    ax.legend()


# -



