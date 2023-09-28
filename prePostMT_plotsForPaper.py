from stable_mt_utils import get_h5
#from matplotlib import rcParams
#rcParams['font.family'] = 'serif'

# +
allmods =! ls ../data/masterOutput//uniformIC/ | grep Default | grep Comp
#print(mods)
print(allmods)

#allmods = allmods[3:] # only look at gammas

#def get_h5(mod):
    
#    return h5.File('../data/betas_and_gammas/uniformIC/{}/COMPAS_Output.h5'.format(mod), 'r')
    
allmods = ['uniformIC/zetaDefault_betaComp_fgamma{}'.format(num) for num in ['0.00', '0.25', '0.50', '0.75', '1.00']]


# +
#myf = h5.File(get_fname(mods[0]))
#myf = get_h5(allmods[4])
#myf.keys()
#MTs = myf['BSE_RLOF']

# +
#mrg = MTs['Merger'][()] == 1
#aPost = MTs['SemiMajorAxis>MT'][()]
#np.sort(aPost[mrg])
#printCompasDetails(MTs, mask=mrg)
# # Decide which parameters to plot
#
# ### For 1st MT
# #### separation, primary stype, secondary stype, primary mass, secondary mass, rad primary, rad secondary 
#

# +
# For objects which enter into stable MT, want to get the rates of mergers, and the separation distribution afterward
# Do this for just the 1st MT episode...
# -

def plotPostStableMtParams(ax, mod, add_title=False, reduceBy=10, bins=20, category=None, color='red', param=0, usePreMT=False): #
    # reduce by means only plot every Nth point, to reduce clutter
    
    myf = get_h5(mod) # TODO add both
            
    #SPs = myf['BSE_System_Parameters']
    MTs = myf['BSE_RLOF']
    #CEs = myf['BSE_Common_Envelopes']
    
    if category not in ['Merger', 'Stable Merger', 'Unstable Merger', 'Stable']:
        print(category)
        raise Exception("no good")
    
    # Get parameters
    mtSeeds = MTs['SEED'][()]
    isCee = MTs['CEE>MT'][()] == 1
    isMerger = MTs['Merger'][()] == 1
    st1 = MTs['Stellar_Type(1)<MT'][()]
    st2 = MTs['Stellar_Type(2)<MT'][()]
    isRlof1 = MTs['RLOF(1)>MT'][()] == 1
    isRlof2 = MTs['RLOF(2)>MT'][()] == 1
    isRlof1Prev = MTs['RLOF(1)<MT'][()] == 1
    isRlof2Prev = MTs['RLOF(2)<MT'][()] == 1
    mask1stMT = np.append(np.array([True]), np.diff(mtSeeds) != 0) # Prepend a true since the first value is always a first seed
    #mask1stMT &= maskPreSN   
    
    ### Extract the primary stellar type at 1st MT if it occurs pre SN 
    mask1stMT = np.append(np.array([True]), np.diff(mtSeeds) != 0) # Prepend a true since the first value is always a first seed
    maskMtContinuation = (isRlof1 & isRlof1Prev) | (isRlof2 & isRlof2Prev) # Need to factor in extended MT
    idxMtContinuation = np.where(maskMtContinuation)[0]
    idxMtContinuationMinus1 = idxMtContinuation - 1
    idxMtStartOfContinuation = idxMtContinuationMinus1[~np.in1d(idxMtContinuationMinus1, idxMtContinuation)]
    idxMtEndOfContinuation = idxMtContinuation[~np.in1d(idxMtContinuation, idxMtContinuationMinus1)]
    maskMtStartOfContinuation = np.in1d(np.arange(len(mtSeeds)),idxMtStartOfContinuation)
    maskMtEndOfContinuation = np.in1d(np.arange(len(mtSeeds)),idxMtEndOfContinuation)
    
    # Ensure everything is preSN
    maskPreSN = ((st1 < 13) | (st1 == 16)) & ((st2 < 13) | (st2 == 16))  # 16 for CHE which are MS stars for this purpose
        
    mask1stMtFirstTimestep = mask1stMT & maskPreSN
    mask1stMtLastTimestep  = (mask1stMtFirstTimestep & ~maskMtStartOfContinuation) 
    # need to add into it, the maskMtend ofcontinuation but only for the ones that apply for the 1st MT
    mask1stMtLastTimestep[maskMtEndOfContinuation] = mask1stMtFirstTimestep[maskMtStartOfContinuation]    
    
    # maskMtStartmerger[maskMtEndOf]
    maskStable1stMt       = np.full(mtSeeds.shape, False)
    mask1stMerger         = np.full(mtSeeds.shape, False)
    mask1stMergerStable   = np.full(mtSeeds.shape, False)
    mask1stMergerUnStable = np.full(mtSeeds.shape, False)

    if usePreMT:
        mask1stMtFirstOrLastTimestep = mask1stMtFirstTimestep
    else:
        mask1stMtFirstOrLastTimestep = mask1stMtLastTimestep
        
    maskStable1stMt[mask1stMtFirstOrLastTimestep]       = (~isCee & ~isMerger)[mask1stMtLastTimestep] #& maskPreSN[mask1stMtFirstTimestep]
    mask1stMerger[mask1stMtFirstOrLastTimestep]         = (isMerger          )[mask1stMtLastTimestep] #& maskPreSN[mask1stMtFirstTimestep]
    mask1stMergerStable[mask1stMtFirstOrLastTimestep]   = (isMerger & ~isCee )[mask1stMtLastTimestep] #& maskPreSN[mask1stMtFirstTimestep]
    mask1stMergerUnStable[mask1stMtFirstOrLastTimestep] = (isMerger & isCee  )[mask1stMtLastTimestep] #& maskPreSN[mask1stMtFirstTimestep]

    
    # Old version
    #maskStable1stMt = mask1stMT & ~isCee & ~isMerger
    #mask1stMerger = mask1stMT & isMerger
    #mask1stMergerStable = mask1stMT & isMerger & ~isCee
    #mask1stMergerUnStable = mask1stMT & isMerger & isCee

        
        
        
    if category == 'Merger':
        maskToUse = mask1stMerger
    elif category == 'Stable Merger':
        maskToUse = mask1stMergerStable
    elif category == 'Unstable Merger':
        maskToUse = mask1stMergerUnstable # set here
    elif category == 'Stable':
        maskToUse = maskStable1stMt # set here
    
    # parameters to plot - post
    aPost = MTs['SemiMajorAxis>MT'][()]
    logAPost = np.log10(aPost)
    st1_mt= MTs['Stellar_Type(1)>MT'][()]
    st2_mt= MTs['Stellar_Type(2)>MT'][()]
    m1_mt = MTs['Mass(1)>MT'][()]
    m2_mt = MTs['Mass(2)>MT'][()]
    r1_mt = MTs['Radius(1)>MT'][()]
    r2_mt = MTs['Radius(2)>MT'][()]
    logR1Post = np.log10(r1_mt)
    logR2Post = np.log10(r2_mt)
    qPost = m2_mt/m1_mt
    
    # pre versions 
    aPre = MTs['SemiMajorAxis<MT'][()]
    logAPre = np.log10(aPre)
    st1_mt_Pre= MTs['Stellar_Type(1)<MT'][()]
    st2_mt_Pre= MTs['Stellar_Type(2)<MT'][()]
    m1_mt_Pre = MTs['Mass(1)<MT'][()]
    m2_mt_Pre = MTs['Mass(2)<MT'][()]
    r1_mt_Pre = MTs['Radius(1)<MT'][()]
    r2_mt_Pre = MTs['Radius(2)<MT'][()]
    logR1Pre = np.log10(r1_mt_Pre)
    logR2Pre = np.log10(r2_mt_Pre)
    qPre = m2_mt_Pre/m1_mt_Pre
    
    if usePreMT:
        paramToPlot = [logAPre , st1_mt_Pre, st2_mt_Pre, m1_mt_Pre , m2_mt_Pre , logR1Pre, logR2Pre, qPre ][param]
    else:
        paramToPlot = [logAPost, st1_mt, st2_mt, m1_mt, m2_mt, logR1Post, logR2Post, qPost][param]
 
    #gamma = float(mod[-3:])
    #label = mod # r'$f_\gamma={:.2f}$'.format(gamma)
    gamma = float(mod[-4:])
    label = r'$f_\gamma={:.2f}$'.format(gamma)

    return paramToPlot[maskToUse], label
    #ax.hist(paramToPlot[maskStable1stMt], bins=bins, histtype='step', label=label, lw=5, color=color)
    #return maskStable1stMt



# USER SETS THIS - choose from the list below


# +
def plotBeforeAndAfter(axs, mods, zeta, param = 0, fs_title = 50, fs_label = 46, fs_corner_letter=20, addLegend=True, corner_letter=None):
    colors = ['red', 'orange', 'green', 'blue', 'violet'] #, 'k'] # for the models
    #mods = allmods[:5]

    ### Top plot should be pre-MT, bottom plot should be (inverted) post-MT

    pltname =  ["logA", "stellarType1", "stellarType2", "mass1", "mass2", "radius1", "radius2", "massRatio"][param]
    title =  ["a pre-/post-MT", "Stellar Type 1 pre-/post-MT", "Stellar Type 2 pre-/post-MT", "Mass 1 pre-/post-MT", "Mass 2 pre-/post-MT", "Radius 1 pre-/post-MT", "Radius 2 pre-/post-MT", "Mass Ratio pre-/post-MT"][param]
    xlabel =  ["$\log_{10}(a / R_\odot)$", "primary stellar type [Hurley]", "secondary stellar type [Hurley]", "primary mass, $M_1 [M_\odot)]$", "secondary mass, $ M_2 [M_\odot)]$", "primary radius, $R_1 [log(R_\odot)]$", "secondary radius, $R_2 [log(R_\odot)]$", "mass ratio $q$" ][param]
    histtype = ['step', 'bar', 'bar', 'step', 'step', 'step', 'step', 'step'][param]

    bins = [ np.linspace(-1.5, 4.5, 25), np.arange(10), np.arange(8), np.linspace(-5, 60, 14), np.linspace(-5, 60, 14), 30, 30, np.linspace(0, 1.2, 25) ][param]

    category = ['Merger', 'Stable Merger', 'Unstable Merger', 'Stable'][3] # choose here

    #mods = allmods
    # Top plot
    ax = axs[0]
    outputs = []
    labels = []
    for ii, mod in enumerate(mods):
        output, label = plotPostStableMtParams(ax, mod=mod, color=colors[ii], param=param, category=category, bins=bins, usePreMT=True)
        outputs.append(output)
        labels.append(label)
    ax.hist(outputs, bins=bins, histtype=histtype, label=labels, lw=5, color=colors)
    
    # Set legend
    if addLegend:
        ax.legend(fontsize=16, title='SMT', title_fontsize=20, loc='upper left')
    
    

    # Set title
    zeta_title = r'$\zeta_{{\mathrm{{{}}}}}$'.format(zeta)
    #ax.set_title("{} for first stable MT episodes :: {}".format(title, zeta_title), fontsize=fs_title)
    #ax.set_title("{}".format(zeta_title), fontsize=fs_title, pad=30)
    ax.set_ylabel("Counts", fontsize=fs_label, y=0)
    ax.tick_params(labelsize=30)

    # Turn ticks off
    #ax.set_xticklabels([])
    ax.label_outer()


    
    
    
    ### Bottom plot
    
    ax = axs[1]
    outputs = []
    labels = []
    for ii, mod in enumerate(mods):
        output, label = plotPostStableMtParams(ax, mod=mod, color=colors[ii], param=param, category=category, bins=bins, usePreMT=False)
        outputs.append(output)
        labels.append(label)
    ax.hist(outputs, bins=bins, histtype=histtype, label=labels, lw=5, color=colors)
    #ax.set_ylabel("Counts", fontsize=fs_label)
    ax.set_xlabel(xlabel, fontsize=fs_label)
    ax.tick_params(labelsize=30)


    if param==0:
        ax.set_xlim(0, 4.5)
    ylimsTop = axs[0].get_ylim()
    ylimsBot = axs[1].get_ylim()

    ylims = np.maximum(ylimsTop, ylimsBot)
    axs[0].set_ylim(*ylims)
    axs[1].set_ylim(*ylims)
    ax.invert_yaxis()


    # Add text explaining pre and post
    xmin, xmax = ax.get_xlim()
    fs_prepost = 40
    axs[0].text(x=xmax-.85, y=ylims[1]*6.4/8, s='Pre-MT', fontsize=fs_prepost)
    axs[1].text(x=xmax-.92, y=ylims[1]*7.4/8, s='Post-MT', fontsize=fs_prepost)
    # Add corner letter
    if corner_letter is not None:
        axs[0].text(x=xmin+.03, y=ylims[1]*6.5/8, s="{})".format(corner_letter), fontsize=fs_corner_letter)
    
if True:
    for ii in range(2):
        #ii = 0 # choose which plot to make

        fig, axes = plt.subplots(nrows=2, figsize=(12, 8), sharex=True)

        axes = axes.flatten()
        zetas = ['Default', 'Default']
        allmods = [['uniformIC/zeta{}_betaComp_fgamma{}'.format(zeta, num) for num in ['0.00', '0.25', '0.50', '0.75', '1.00']] for zeta in zetas]

        plotBeforeAndAfter(axes, allmods[ii], zeta=zetas[ii], addLegend=False)

        plt.rcParams.update({"savefig.facecolor": 'white', })
        fig.tight_layout()
        plt.savefig('../plots/final_plots/aPreAPostMirror_zeta{}.pdf'.format(zetas[ii]), format='pdf', bbox_inches='tight')


# +
def plotPostStableMtSepRatio(ax, mods, zeta, param = 0,  fs_title = 50, fs_label = 40, addYlabel=True, fs_corner_letter=20, corner_letter=None):
    colors = ['red', 'orange', 'green', 'blue', 'violet'] #, 'k'] # for the models
    #mods = allmods[:5]

    ### Top plot should be pre-MT, bottom plot should be (inverted) post-MT

    pltname =  ["logA", "stellarType1", "stellarType2", "mass1", "mass2", "radius1", "radius2", "massRatio"][param]
    title =  ["a pre-/post-MT", "Stellar Type 1 pre-/post-MT", "Stellar Type 2 pre-/post-MT", "Mass 1 pre-/post-MT", "Mass 2 pre-/post-MT", "Radius 1 pre-/post-MT", "Radius 2 pre-/post-MT", "Mass Ratio pre-/post-MT"][param]
    xlabel =  ["Semi-major Axis, $\log(a / R_\odot)$", "primary stellar type [Hurley]", "secondary stellar type [Hurley]", "primary mass, $M_1 [M_\odot)]$", "secondary mass, $ M_2 [M_\odot)]$", "primary radius, $R_1 [log(R_\odot)]$", "secondary radius, $R_2 [log(R_\odot)]$", "mass ratio $q$" ][param]
    histtype = ['step', 'bar', 'bar', 'step', 'step', 'step', 'step', 'step'][param]

    bins = [ np.linspace(0, 3, 25), np.arange(10), np.arange(8), np.linspace(-5, 60, 14), np.linspace(-5, 60, 14), 30, 30, np.linspace(0, 1.2, 25) ][param]

    category = ['Merger', 'Stable Merger', 'Unstable Merger', 'Stable'][3] # choose here

    #mods = allmods
    # Top plot
    #ax = axs[0]
    outputs = []
    labels = []
    for ii, mod in enumerate(mods):
        logAPre, label = plotPostStableMtParams(ax, mod=mod, color=colors[ii], param=param, category=category, bins=bins, usePreMT=True)
        logAPost, label = plotPostStableMtParams(ax, mod=mod, color=colors[ii], param=param, category=category, bins=bins, usePreMT=False)

        #aRatio = np.power(10, logAPost)/np.power(10, logAPre)
        aRatio = np.power(10, logAPost)/np.power(10, logAPre)
        outputs.append(aRatio)
        labels.append(label)
    ax.hist(outputs, bins=bins, histtype=histtype, label=labels, lw=5, color=colors)
    
    # Set legend
    ax.legend(fontsize=24, title='SMT', title_fontsize=30)
    

    # Set title
    zeta_title = r'$\zeta_{{\mathrm{{{}}}}}$'.format(zeta)
    #ax.set_title("{} for first stable MT episodes :: {}".format(title, zeta_title), fontsize=fs_title)
    #ax.set_title("{}".format(zeta_title), fontsize=fs_title, pad=30)
    if addYlabel:
        ax.set_ylabel("Counts", fontsize=fs_label)
    ax.set_xlabel("Orbital tightening $a_{\mathrm{post}}/a_{\mathrm{pre}}$", fontsize=fs_label)
    #ax.set_ylim(0, 9000)
    ax.set_xlim(-.3, 3.1)
    
    # Add corner letter
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if corner_letter is not None:
        ax.text(x=xmin+.03, y=ymax*7.2/8, s="{})".format(corner_letter), fontsize=fs_corner_letter)

    
#ii = 0
#for ii in range(2):
if False:
    fig, ax = plt.subplots(figsize=(8, 6))#, sharex=True)

    #axes = axes.flatten()
    zetas = ['Default', 'Default']
    allmods = [['uniformIC/zeta{}_betaComp_fgamma{}'.format(zeta, num) for num in ['0.00', '0.25', '0.50', '0.75', '1.00']] for zeta in zetas]
    addY = [True, True]

    #for ii, ax in enumerate(axes):
    #ax = axes
    plotPostStableMtSepRatio(ax, allmods[ii], zeta=zetas[ii], addYlabel=addY[ii])

    plt.rcParams.update({"savefig.facecolor": 'white', })
    fig.tight_layout()
    plt.savefig('../plots/final_plots/aPreDividedByAPost_zeta{}.pdf'.format(zetas[ii]), format='pdf', bbox_inches='tight')
# +
def plotPostStableMtSepRatioLog(ax, mods, zeta, param = 0,  fs_title = 50, fs_label = 46, addYlabel=True, fs_corner_letter=20, corner_letter=None):
    colors = ['red', 'orange', 'green', 'blue', 'violet'] #, 'k'] # for the models
    #mods = allmods[:5]

    ### Top plot should be pre-MT, bottom plot should be (inverted) post-MT

    pltname =  ["logA", "stellarType1", "stellarType2", "mass1", "mass2", "radius1", "radius2", "massRatio"][param]
    title =  ["a pre-/post-MT", "Stellar Type 1 pre-/post-MT", "Stellar Type 2 pre-/post-MT", "Mass 1 pre-/post-MT", "Mass 2 pre-/post-MT", "Radius 1 pre-/post-MT", "Radius 2 pre-/post-MT", "Mass Ratio pre-/post-MT"][param]
    xlabel =  ["Semi-major Axis, $\log(a / R_\odot)$", "primary stellar type [Hurley]", "secondary stellar type [Hurley]", "primary mass, $M_1 [M_\odot)]$", "secondary mass, $ M_2 [M_\odot)]$", "primary radius, $R_1 [log(R_\odot)]$", "secondary radius, $R_2 [log(R_\odot)]$", "mass ratio $q$" ][param]
    histtype = ['step', 'bar', 'bar', 'step', 'step', 'step', 'step', 'step'][param]

    #bins = [ np.linspace(0, 3, 25), np.arange(10), np.arange(8), np.linspace(-5, 60, 14), np.linspace(-5, 60, 14), 30, 30, np.linspace(0, 1.2, 25) ][param]
    bins = [ np.linspace(-2.5, 1.25, 31), np.arange(10), np.arange(8), np.linspace(-5, 60, 14), np.linspace(-5, 60, 14), 30, 30, np.linspace(0, 1.2, 25) ][param]

    category = ['Merger', 'Stable Merger', 'Unstable Merger', 'Stable'][3] # choose here

    #mods = allmods
    # Top plot
    #ax = axs[0]
    outputs = []
    labels = []
    for ii, mod in enumerate(mods):
        logAPre, label = plotPostStableMtParams(ax, mod=mod, color=colors[ii], param=param, category=category, bins=bins, usePreMT=True)
        logAPost, label = plotPostStableMtParams(ax, mod=mod, color=colors[ii], param=param, category=category, bins=bins, usePreMT=False)

        #aRatio = np.power(10, logAPost)/np.power(10, logAPre)
        aRatio = np.power(10, logAPost)/np.power(10, logAPre)
        logARatio = np.log10(aRatio)
        print(np.sort(logARatio)[0])
        print(np.sort(logARatio)[-1])
        outputs.append(logARatio)
        labels.append(label)
    
    
    
    ax.hist(outputs, bins=bins, histtype=histtype, label=labels, lw=5, color=colors)
    
    # Set legend
    ax.legend(fontsize=32, title='AM loss', title_fontsize=40)
    

    # Set title
    zeta_title = r'$\zeta_{{\mathrm{{{}}}}}$'.format(zeta)
    #ax.set_title("{} for first stable MT episodes :: {}".format(title, zeta_title), fontsize=fs_title)
    #ax.set_title("{}".format(zeta_title), fontsize=fs_title, pad=30)
    if addYlabel:
        ax.set_ylabel("Counts", fontsize=fs_label)
    ax.set_xlabel("$\log_{10}(a_{\mathrm{post}}/a_{\mathrm{pre}})$", fontsize=fs_label)
    #ax.set_ylim(0, 9000)
    ax.set_xlim(-2.5, 1.2)
    ax.tick_params(labelsize=30)

    
    # Add corner letter
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if corner_letter is not None:
        ax.text(x=xmin+.03, y=ymax*7.2/8, s="{})".format(corner_letter), fontsize=fs_corner_letter)

    
#ii = 0
#for ii in range(2):
if True:
    fig, ax = plt.subplots(figsize=(12, 8))#, sharex=True)

    #axes = axes.flatten()
    zetas = ['Default', 'Default']
    allmods = [['uniformIC/zeta{}_betaComp_fgamma{}'.format(zeta, num) for num in ['0.00', '0.25', '0.50', '0.75', '1.00']] for zeta in zetas]
    addY = [True, True]

    #for ii, ax in enumerate(axes):
    #ax = axes
    plotPostStableMtSepRatioLog(ax, allmods[ii], zeta=zetas[ii], addYlabel=addY[ii])

    plt.rcParams.update({"savefig.facecolor": 'white', })
    fig.tight_layout()
    plt.savefig('../plots/final_plots/aPreDividedByAPost_zeta_log{}.pdf'.format(zetas[ii]), format='pdf', bbox_inches='tight')
# -




# +
# Plot the mirrored pre/post and the ratio post/pre in a 3 plot column (for a given zeta)

# Mirror plot
#ii = 0 # choose which plot to make
#fig, axes = plt.subplots(nrows=3, figsize=(12, 8))#, sharex=True)
for ii in range(2):
    moddir = ['uniformIC', 'moedistef'][ii] #'uniformIC', , 'moedistef'
    ii = ii%2
    allmods = ['{}/zetaDefault_betaComp_fgamma{}'.format(moddir, num) for num in ['0.00', '0.25', '0.50', '0.75', '1.00']]# for zeta in zetas]
    print(allmods)

    fig = plt.figure(constrained_layout=False, figsize=(12, 15))

    gs = fig.add_gridspec(nrows=2, ncols=1, hspace=.3)
    gs0 = gs[0].subgridspec(2, 1, hspace=0) 
    gs1 = gs[1] # single

    ax1 = fig.add_subplot(gs0[0])
    ax2 = fig.add_subplot(gs0[1], sharex=ax1)
    ax3 = fig.add_subplot(gs1)

    #corner_letters = [['i', 'iii'],['ii', 'iv']][ii]
    corner_letters = [None, None]
    fs_corner_letter = 24

    axs = [ax1, ax2]
    zetas = ['Default']#, 'Ge20']

    plotBeforeAndAfter(axs, allmods, zeta=zetas[0], addLegend=False, fs_corner_letter=fs_corner_letter)

    addY = [True, True]

    #ax = axes[2]
    plotPostStableMtSepRatio(ax3, allmods, zeta=zetas[0], addYlabel=addY[ii], fs_corner_letter=fs_corner_letter)

    plt.rcParams.update({"savefig.facecolor": 'white', })
    #fig.tight_layout()
    plt.savefig('../plots/final_plots/aPreAPostTriplet_{}_zeta{}.pdf'.format(moddir, zetas[0]), format='pdf', bbox_inches='tight')
# -










