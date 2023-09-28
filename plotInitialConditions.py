# # Plot initial parameters as scatter, q vs a, with size given by m1, 
# # colored by end state of binary

# ### TODO: add in CEE and merger as same endpoint, and 1st SN as another



# +
# TODO: End points are stable MT, unstable MT, and non-interactive (by 1st SN)

# Need to get some way to identify first MT, and relative to first SN. 

# Color based on stellar type of 1st interaction

# triangle for cee/merger (before 1st SN), circle for non-merger (before 1st SN)

# Need to create a mask for first MT - then check stellar types (of both) to make sure it's preSN, as well as the primary stype for color
# Check whether it merged at all preSN (in SPs)

# +
# allmods =!ls ../data/
allmods = allmods[3:] # only look at gammas
print(allmods)

def get_h5(mod):
    return h5.File('../data/{}/COMPAS_Output.h5'.format(mod), 'r')
    


# -

myf = get_h5(allmods[1])
SPs = myf['BSE_System_Parameters']
MTs = myf['BSE_RLOF']
printCompasDetails(MTs)


def plotInitialConditionDistribution(ax, mod, add_title=False, title=None, reduceBy=10):
    # reduce by means only plot every Nth point, to reduce clutter
    

    myf = get_h5(mod) # TODO add both

    #myf = get_h5(mod) # TODO add both
    SPs = myf['BSE_System_Parameters']
    MTs = myf['BSE_RLOF']
    CEs = myf['BSE_Common_Envelopes']
    
    # Get parameters
    spSeeds = SPs['SEED'][()]
    mtSeeds = MTs['SEED'][()]
    ceSeeds = CEs['SEED'][()]
    m1 = SPs['Mass@ZAMS(1)'][()]
    m2 = SPs['Mass@ZAMS(2)'][()]
    st1_f = SPs['Stellar_Type(1)'][()]
    st2_f = SPs['Stellar_Type(2)'][()]
    st1_mt= MTs['Stellar_Type(1)<MT'][()]
    st2_mt= MTs['Stellar_Type(2)<MT'][()]
    aZams = SPs['SemiMajorAxis@ZAMS'][()] # AU
    eZams = SPs['Eccentricity@ZAMS'][()]    
    st1_ce = CEs['Stellar_Type(1)<CE'][()]
    st2_ce = CEs['Stellar_Type(2)<CE'][()]
    
    # Derived quantities
    q = m2/m1
    rP = aZams*(1-eZams)
    
    # Define resizing functions for the size of the dots
    def resizeM(mvals):
        mvs = mvals.copy()
        return mvs*mvs/5
    mm1 = resizeM(m1)


    ### Color based on the primary stellar type at 1st MT if it occurs pre SN (need a color/mask for non-interacting)
    mask1stMT = np.append(np.array([True]), np.diff(mtSeeds) != 0) # Prepend a true since the first value is always a first seed
    maskPreSN = (st1_mt < 13) & (st2_mt < 13)
    mask1stMtPreSN = mask1stMT & maskPreSN
    st1mt_1stPreSnMt = st1_mt[mask1stMtPreSN]
    # Need to remap these types back onto a SP-like array
    sdsMt_1stPreSnMt = mtSeeds[mask1stMtPreSN]
    spMask_1stPreSnMt = np.in1d(spSeeds, sdsMt_1stPreSnMt)
    st1sp_1stPreSnMt = np.zeros_like(spSeeds)                # initialize all stypes to 0, corresponding to no interactions
    st1sp_1stPreSnMt[spMask_1stPreSnMt] = st1mt_1stPreSnMt   # update the interactors to the correct stype
    #print(np.unique(stype1_1stPreSnMt, return_counts=True)) # 1-6
    
    

    ### Shape based on CEE/Merger, interacting non-merger, and isolated

    # Get Merger data from SPs
    isMerger = SPs['Merger'][()] == 1
    ends_preSN = (st1_f < 13) & (st2_f < 13)
    maskMergesPreSn = isMerger & ends_preSN # don't include post SN mergers
    # Get CEE data from CEs
    maskCePreSN = (st1_ce < 13) & (st2_ce < 13)
    sdsCePreSN = ceSeeds[maskCePreSN]
    maskSp_CePreSN = np.in1d(spSeeds, sdsCePreSN)
    # Combine mergers and CEE
    maskMergerOrCeePreSn = maskMergesPreSn | maskSp_CePreSN
    
    # Get interacting non-mergers by grabbing all interacters preSN, and removing any from the merger mask above
    maskMtInteractionsPreSn = (st1_mt < 13) & (st2_mt < 13) # all MT events that occur prior to a SN
    sdsMtInteractsPreSn = np.unique(mtSeeds[maskMtInteractionsPreSn]) # the seeds of everything that interacts prior to a SN
    maskSpInteractPreSn = np.in1d(spSeeds, sdsMtInteractsPreSn)
    maskSpInteractNonMerger = maskSpInteractPreSn & ~maskMergerOrCeePreSn
    
    # Everything else is a non-interactor - combine the masks together
    state_masks = [maskMergerOrCeePreSn, maskSpInteractNonMerger, ~maskMergerOrCeePreSn & ~maskSpInteractNonMerger]
    # check the masks
    for mm in range(3):
        msk1 = state_masks[mm%3]
        msk2 = state_masks[(mm+1)%3]
        assert ~np.any(msk1 & msk2)
    
    # Shape based on CEE/Merger, interacting, and isolated
    shapes = ['*', '^', 'o']
    shape_labels = ['Merger | CEE', 'Interacting', 'Isolated']
    # Color based on stellar type of 1st interaction - black for no interaction - color_labels are just the stellar type integers
    colors = ['black', 'darkorange', 'red', 'blue', 'green', 'sandybrown', 'fuchsia']
    
    
    
    ### Additional masks to help plotting
    reduction = (spSeeds%reduceBy == 0) # reduce number of points for clarity - leave this
    massBounds = (m1 < 40) & (m1 > 10) # don't care about very massive stars or very low mass stars - adjust this as needed
    
    # Plot the points
    for ii, state_mask in enumerate(state_masks):

        state_mask &= reduction & massBounds
        shape = shapes[ii]

        for stype in range(7):
            #print(stype) # check it goes 0 through 6
            
            # Apply additional mask on stype, for the color
            stype_mask = st1sp_1stPreSnMt == stype
            color = colors[stype]
            
            mask = state_mask & stype_mask
            ax.scatter(rP[mask], q[mask], mm1[mask], marker=shape, c=color)
            
            if False: # show seed values for these points
                strictmask = (q>.85)&(q<.95)&(rP>.2)&(rP<.4)
                mask2 = strictmask&mask
                for rpp, qqq, spsds in zip(rP[mask2], q[mask2], spSeeds[mask2]):
                    ax.text(rpp, qqq-.05, spsds)
               
            
        
        
        
    # Set plot attributes
    ax.set_xlabel('Periapsis (AU)', fontsize=18)
    ax.set_xscale('log')
    ax.set_xlim(.01, int(4e4)) 
    ax.set_ylabel(r'Mass Ratio $q$', fontsize=18)
    ax.set_yscale('linear')
    ax.tick_params(labelsize=18)

    if add_title:
        if title is None:            
            gamma = float(mod[-3:])
            title = r'$f_\gamma={:.2f}$'.format(gamma)
        ax.set_title(title, fontsize=24)
        #ax.set_title(mod, fontsize=18)
    
    
    
    
    ### Produce custom legends showing shapes, sizes, and colors individually
    
    # Add custom legend for shapes
    custom_legend = []
    for ii in range(3):
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=shape_labels[ii], color='k', alpha=1.0, lw=0, marker=shapes[ii], ms=18))
    lgnd1 = ax.legend(handles=custom_legend, title='End state (pre-SN)', fontsize=14, title_fontsize=18, loc='upper right')
    ax.add_artist(lgnd1)

    
    # Add custom legend for colors
    custom_legend = []
    for stype in range(7):
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=stype, color=colors[stype], alpha=1.0, lw=0, marker='s', ms=18))
    lgnd2 = ax.legend(handles=custom_legend, title='Primary type at 1st interaction', fontsize=14, title_fontsize=18, loc='lower right', ncol=2)
    ax.add_artist(lgnd2)

    
    #for col, lab in zip(colors, labels):      
    #    custom_legend.append(mpl.lines.Line2D([0], [0], color=col, label=lab, lw=0, ms=10, marker='o'))
    
    # Add custom legend for masses
    custom_handles = []
    custom_labels = []
    for m_str in ['10', '25', '40']:
        m_val = float(m_str)
        resizedM=np.sqrt(resizeM(np.array([m_val])))
        l1 = mpl.lines.Line2D([0, 0], [0, 0], ms=resizedM, color='k', alpha=0.5, lw=0, marker=shapes[0])
        l2 = mpl.lines.Line2D([0, 10], [0, 0], ms=resizedM, color='k', alpha=0.5, lw=0, marker=shapes[1])
        l3 = mpl.lines.Line2D([50], [0], ms=resizedM, color='k', alpha=0.5, lw=0, marker=shapes[2])
        custom_handles.append((l1, l2, l3))
        custom_labels.append(r"$M_1={}M_\odot$".format(m_str))
    lgnd2 = ax.legend(handles=custom_handles, labels=custom_labels, handler_map={tuple: mpl.legend_handler.HandlerTuple(None)}, 
              title='Primary mass at ZAMS', fontsize=10, title_fontsize=18, prop={'size':25}, loc='center right')



    # Add demarkations for the various sub-models
    if False: #(not plotm1vsRp) and (mod == 'fiducial'):
        ax.plot([.02, 4e3], [.5, .5], 'k')
        ax.plot([10, 10], [.05, 1.08], 'k')
        ax.text(x=.3, y=1.03, s='aMax10AU', fontsize=fs_minmaxmods)
        ax.text(x=30, y=1.03, s='aMin10AU', fontsize=fs_minmaxmods)
        ax.text(x=1.1e3, y=0.75, s='qMin0.5', fontsize=fs_minmaxmods)
        ax.text(x=1.1e3, y=0.25, s='qMax0.5', fontsize=fs_minmaxmods)


# +
mods = [allmods[0]] + [allmods[2]] + [allmods[-1]]
print(mods)
fig, axes = plt.subplots(nrows=len(mods), figsize=(25, 28))

for ii, mod in enumerate(mods):
    ax = axes[ii]
    plotInitialConditionDistribution(ax, mod, add_title=True)
    

#plt.rcParams.update({
#"savefig.facecolor": 'white', # (0.0, 0.0, 1.0, 0.2),  # blue  with alpha = 20%
#})

#plt.savefig('../data/plots/endstatesFromICs_g0vsg1.png')


# +
newmods =! ls ../data/new_runs/moedistef/
moddir = '../data/new_runs/moedistef/'
print(newmods)

fig, axes = plt.subplots(nrows=len(newmods), figsize=(25, 100))

for ii, mod in enumerate(newmods):
    ax = axes.flatten()[ii]
    plotInitialConditionDistribution(ax, moddir+mod, add_title=True, title=mod, reduceBy=10)


# -






