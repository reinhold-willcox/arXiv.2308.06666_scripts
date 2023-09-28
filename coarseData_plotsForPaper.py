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

# # This file contains all the individual plots to process the data as I see fit. 
#
# ## These might be combined or resized as needed later, so they should be built flexibly

import numpy as np
import matplotlib as mpl
from stable_mt_utils_fix import get_h5, extractEndStates

# +
colorBlindFriendyPalette = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

fig, ax = plt.subplots()

for ii in range(len(colorBlindFriendyPalette)):
    ax.plot([0, 1], [ii, ii], color=colorBlindFriendyPalette[ii], lw=8)
# -

# testing
testing = False
if testing:
    moddirs = ['moedistef', 'uniformIC']
    mod = 'single_variations/' + moddirs[0] + '/beta0.00/'




# +
def makePQplot(ax, mod=None, Data=None, add_title=False, title=None, reduceBy=10, includeSingles=True, useEndState=True, fs_axis_label=24, fs_title=30, lg_title=30, lg_fs=24, leg_lw=2, xUpp=None,  **kwargs):
    
    if (mod == None) == (Data == None):
        raise Exception("Need to submit exactly one of mod or Data")
        
    ### Extract the post-processed data
    params, stype_masks_colors, state_masks_markers_labels_colors = \
        extractEndStates(mod, **kwargs) if Data is None else Data

    [q, rP, m1, mT] = params
    [stype_masks, stype_colors] = stype_masks_colors
    [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors
    
    # Choose whether to plot the mt1_state or the end_state
    if useEndState:
        state_masks = end_state_masks 
    else:
        state_masks = mt1_state_masks 

    
    ### Additional masks to help plotting
    reduction = (np.arange(len(m1)) % reduceBy == 0) # reduce number of points for clarity - leave this
    massBounds = (m1 < 40) & (m1 > 10) # don't care about very massive stars or very low mass stars - adjust this as needed

    # Define resizing functions for the size of the dots
    f10=100
    f40=1000
    fAlpha=(f40-f10)/30
    fBeta = (4*f10-f40)/3
    def resizeM(mvals):
        mvs = mvals.copy()
        return fAlpha*mvs+fBeta #mvs*mvs/2
    mm1 = resizeM(m1)

        
    ############################
    ### Scatter plot
    for stateId in np.unique(state_masks):
        
        state_mask = (state_masks==stateId) & reduction & massBounds
        marker = state_markers[stateId]

        for stype in np.unique(stype_masks):
            
            if stype > 6:
                continue
            # Apply additional mask on stype, for the color
            stype_mask = (stype_masks==stype)
            color = stype_colors[stype]

            mask = state_mask & stype_mask & (q>0.1)
            ax.scatter(rP[mask], q[mask], mm1[mask], marker=marker, c=color)

            if False: # show seed values for these points - for debugging
                strictmask = (q>.85)&(q<.95)&(rP>.2)&(rP<.4)
                mask2 = strictmask & mask
                for rpp, qqq, spsds in zip(rP[mask2], q[mask2], spSeeds[mask2]):
                    ax.text(rpp, qqq-.05, spsds)

    # Set plot attributes
    ax.set_xlabel(r'Periapsis $r_P$ [AU]', fontsize=fs_axis_label)
    ax.set_xscale('log')
    ax.set_ylim(-.05, 1.05)
    if xUpp is None:
        xUpp = 3e2 if includeSingles else 1e2
    ax.set_xlim(.03, int(xUpp)) 
    ax.set_ylabel(r'Mass Ratio $q$', fontsize=fs_axis_label)
    ax.set_yscale('linear')
    ax.tick_params(labelsize=fs_axis_label-10)


    ### Produce custom legends showing markers, sizes, and colors individually

    
    # Add custom legend for markers
    custom_legend = []
    lgnd_title = 'End state' if useEndState else 'Outcome of 1st MT'
    for stateId in np.unique(state_masks):
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=state_labels[stateId], color='k', alpha=1.0, lw=0, marker=state_markers[stateId], ms=24))
    lgnd1 = ax.legend(handles=custom_legend, title=lgnd_title, fontsize=lg_fs, title_fontsize=lg_title, framealpha=1, loc=(.67, .43), ncol=2, handletextpad=.01, columnspacing=.8)
    lgnd1.get_frame().set_linewidth(leg_lw)
    ax.add_artist(lgnd1)


    # Add custom legend for colors
    custom_legend = []
    stype_num_name = {0: 'N/A', 1: 'MS', 2: 'HG', 3:'GB', 4:'CHeB', 5:'AGB'}#, 6:'TPAGB'}
    for stype in np.arange(6): # np.unique(stype_masks)[:-1]:
        if (stype == 0) and not (includeSingles):
            continue
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=stype_num_name[stype], color=stype_colors[stype], alpha=1.0, lw=0, marker='s', ms=24))
    lgnd2 = ax.legend(handles=custom_legend, title='Primary type at 1st MT', fontsize=lg_fs, title_fontsize=lg_title, framealpha=1, loc='upper right', ncol=2)
    lgnd2.get_frame().set_linewidth(leg_lw)
    ax.add_artist(lgnd2)
    

    # Add custom legend for masses
    custom_handles = []
    custom_labels = []
    if includeSingles:
        XX_YY = [([0, 0], [0, 0]), 
                 ([0, 10], [0, 0]),
                 ([50], [0]),
                 ([20], [0])]           
    else:
        XX_YY = [([0, 0], [0, 0]), 
                 ([0, 10], [0, 0]),
                 ([50], [0])]           
    for m_str in ['10', '25', '40']:
        m_val = float(m_str)
        resizedM=np.sqrt(resizeM(np.array([m_val])))
        custom_tuple = []
        for jj, xx_yy in enumerate(XX_YY):
            custom_tuple.append(mpl.lines.Line2D(xx_yy[0], xx_yy[1], ms=resizedM, color='k', alpha=0.5, lw=0, marker=state_markers[jj]))
        custom_handles.append(tuple(custom_tuple))
        custom_labels.append(r"$M_1={}M_\odot$".format(m_str))
    lgnd3 = ax.legend(handles=custom_handles, labels=custom_labels, handler_map={tuple: mpl.legend_handler.HandlerTuple(None)}, \
              title='Primary mass at ZAMS', fontsize=8, framealpha=1, title_fontsize=lg_title, prop={'size':50}, loc='lower right')
    lgnd3.get_frame().set_linewidth(leg_lw)
    
    if title is not None:
        add_title = True
    if add_title:
        if title is None: 
            modcat, moddir, modname = mod.split('/')
            title = "{} :: {} :: {}".format(modcat, moddir, modname)
        ax.set_title(title, fontsize=fs_title)
    
#testing
testing = False
#'../data/masterOutput/uniformIC/'
mod = 'moedistef/qCritGe20'
if testing:
    fig, ax = plt.subplots(figsize=(25, 8))
    makePQplot(ax, mod, add_title=False, title=None, reduceBy=20, fs_axis_label=24)
# -




# +

def plotUniformVsMoeDiStefPQsOnly(savedir=None):

    # Want to compare P-q plots, with singles, with the pie chart on the side
    
    moddirs = ['uniformIC', 'moedistef']
    mods = [moddir + '/zetaDefault_betaComp_fgamma0.00/' for moddir in moddirs]
    
    fig, axes = plt.subplots(nrows=2, figsize=(25, 30))

    # Get Data
    Data = [extractEndStates(mod=mod) for mod in mods]
    
    xUpp = 1e4
    reduceBy = 30
    titles = [ 'Traditional Sampler', 'Moe \& Di Stefano (2017)' ]
    for ii in range(2):
        makePQplot(ax=axes[ii], Data=Data[ii], title=titles[ii], xUpp=xUpp, reduceBy=reduceBy, fs_axis_label=50, fs_title=60, lg_title=50, lg_fs=40, leg_lw=4, useEndState=False)    
    
    fig.tight_layout()
    
    
# testing
plotUniformVsMoeDiStefPQsOnly()

plt.savefig('../plots/final_plots/{}.pdf'.format('uniformIC_vs_moeDiStef_PQs_only'), format='pdf')

# -






# +
def makeHistogram(ax, mod=None, Data=None, which=None, useEndState=True, add_title=False, title=None, fs_axis_label=24, **kwargs):
            
    if (mod == None) == (Data == None):
        raise Exception("Need to submit exactly one of mod or Data")
        
    ### Extract the post-processed data
    params, stype_masks_colors, state_masks_markers_labels_colors = \
        extractEndStates(mod, **kwargs) if Data is None else Data

    [q, rP, m1, mT] = params
    [stype_masks, stype_colors] = stype_masks_colors
    [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors
     
    # Choose whether to plot the mt1_state or the end_state
    if useEndState:
        state_masks = end_state_masks 
    else:
        state_masks = mt1_state_masks 
    

    which_options = ['q', 'rP', 'm1', 'mT']    
    if which not in which_options:
        raise Exception("Need to choose a param to plot, one of: " + ', '.join(which_options))
        return
    else:
        ii = which_options.index(which)
    
    
    #masks = state_masks
    labels = state_labels[:len(state_masks)]
    labels.append('Total')
    colors = ['red', 'green', 'blue', 'black']
    ls = ['--', ':', '-.', '-']

    
    param =  [q,   np.log10(rP),                                 np.log10(m1),                  np.log10(mT)][ii]
    xlabel = ['Mass Ratio at ZAMS, $q$', r'Periapsis at ZAMS, $\log_{{10}}$($r_P$/AU)', r'Primary Mass at ZAMS, $\log_{{10}}(M_1/M_\odot)$', r'Total Mass at ZAMS, $\log_{{10}}(M_T/M_\odot)$'][ii]   

    data = []
    for stateId in np.unique(state_masks):
        if stateId == 0:
            continue # don't include singles in the histograms
        data.append(param[state_masks == (stateId)])
    data.append(param[state_masks != 0]) # Get all non-singles 
    
    bins = [
        np.linspace(0.05, 1.05, 41),
        np.linspace(-1.5, 1.2, 55),
        np.linspace(1, 2, 41),
        np.linspace(1, 2, 41),
    ][ii]

    for jj in range(len(data)):
        ax.hist(data[jj], bins=bins, color=colors[jj], label=labels[jj], histtype='step', lw=3, linestyle=ls[jj])

    ax.set_xlim((bins[0], bins[-1]))
    ax.set_xlabel(xlabel, fontsize=fs_axis_label)
    ax.set_ylabel('Counts', fontsize=fs_axis_label)

    title = 'End state (pre-SN)' if useEndState else 'First MT outcome'
    ax.legend(title=title, fontsize=18, title_fontsize=24, framealpha=1, loc='upper right')
        
    
        
#testing
if False:
    fig, ax = plt.subplots(figsize=(15, 8))
    #  which_options = ['q', 'rP', 'm1', 'mT']    
    makeHistogram(ax, mod, which='m1', useEndState=False, fs_axis_label=24)        


# +
def makePieChart(ax, mod=None, Data=None, add_title=False, useEndState=True, includeSingles=False, title=None, labeldistance=1.1, **kwargs):
            
    if (mod == None) == (Data == None):
        raise Exception("Need to submit exactly one of mod or Data")
        
    ### Extract the post-processed data
    params, stype_masks_colors, state_masks_markers_labels_colors = \
        extractEndStates(mod, **kwargs) if Data is None else Data

    [q, rP, m1, mT] = params
    [stype_masks, stype_colors] = stype_masks_colors
    [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors
     
    # Choose whether to plot the mt1_state or the end_state
    if useEndState:
        state_masks = end_state_masks 
    else:
        state_masks = mt1_state_masks 
    
    data = []
    for stateId in np.unique(state_masks):
        if (not includeSingles) and (stateId == 0):
            continue # don't include singles
        data.append(np.sum([state_masks == (stateId)]))
    
    
    fracs = data/np.sum(data)
    labels = ["{}\n{:.0f}\%".format(lab, frac*100) for lab, frac in zip(state_labels[:len(state_masks)], fracs)]
    
    colors = state_colors
    if not includeSingles:
        colors = colors[1:]
        
    ax.pie(data, colors=colors, labels=labels, labeldistance=labeldistance, textprops={'fontsize':24}, explode=[.1, 0, 0, .1])    
    
    if add_title:
        ax.set_title(title, fontsize=40, pad=50)
        
#testing
if testing:
    
    fig, axes = plt.subplots(ncols=2, figsize=(15, 8))
    mods = ['uniformIC/fiducial/', 'moedistef/fiducial/']
    titles = ['Traditional Sampler', 'Moe \& Di Stefano (2017)']
    for ii, ax in enumerate(axes):
        makePieChart(ax, mods[ii], includeSingles=True, add_title=True, title=titles[ii], labeldistance=1.2, useEndState=False)
        
    fig.tight_layout()
    plt.savefig('../plots/final_plots/uniform_vs_modistefano_piecharts.pdf', format='pdf', bbox_inches='tight')

    
    

    
# -




# +
def makePieChart2(ax, mod=None, Data=None, add_title=False, useEndState=True, includeSingles=False, title=None, labeldistance=1.1, ii=None, **kwargs):
            
    if (mod == None) == (Data == None):
        raise Exception("Need to submit exactly one of mod or Data")
        
    ### Extract the post-processed data
    params, stype_masks_colors, state_masks_markers_labels_colors = \
        extractEndStates(mod, **kwargs) if Data is None else Data

    [q, rP, m1, mT] = params
    [stype_masks, stype_colors] = stype_masks_colors
    [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors
     
    # Choose whether to plot the mt1_state or the end_state
    if useEndState:
        state_masks = end_state_masks 
    else:
        state_masks = mt1_state_masks 
    
    data = []
    for stypeId in np.unique(stype_masks):
        if (not includeSingles) and (stypeId == 0):
            continue # don't include singles
        data.append(np.sum([(stype_masks == (stypeId)) & (state_masks==ii)]))
        
    print(state_labels[ii])
    
    stype_num_name = {0: 'N/A', 1: 'MS', 2: 'HG', 3:'GB', 4:'CHeB', 5:'AGB', 6:'TPAGB'}
    stype_labels = list(stype_num_name.values())
    
    fracs = data/np.sum(data)
    labels = ["{}\n{:.0f}\%".format(lab, frac*100) for lab, frac in zip(stype_labels[:len(stype_masks)], fracs)]
    
    colors = state_colors
    if not includeSingles:
        colors = colors[1:]
        
    ax.pie(data,  labels=labels)#, labeldistance=labeldistance, textprops={'fontsize':24}, explode=[.1, 0, 0, .1])    
    
    if add_title:
        ax.set_title(title, fontsize=40, pad=50)
        
#testing
if False:
    
    fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(15, 20))
    mods = ['uniformIC/fiducial/', 'moedistef/fiducial/']
    titles = ['Traditional Sampler', 'Moe \& Di Stefano (2017)']
    for ii, ax in enumerate(axes.flatten()):
        makePieChart2(ax, mods[ii%2], includeSingles=True, add_title=bool(ii<2), title=titles[ii%2], labeldistance=1.2, useEndState=False, ii=ii//2+1)
        
    fig.tight_layout()
    #plt.savefig('../plots/final_plots/uniform_vs_modistefano_piecharts.pdf', format='pdf', bbox_inches='tight')

    
    

    
# -










# +
def getTstring(Tconv):
    return r'$T_{{\mathrm{{conv, \;{}}}}}$'.format(Tconv[5:])
def getBetaString(beta):
    bstr = beta[4:]
    if bstr == 'Comp':
        return r'$\beta_{\mathrm{Comp}}$' 
    else:
        return r'$\beta={}$'.format(bstr)
def getFGammaString(fgamma):
    return r'$f_{\gamma}=$' + fgamma[6:]
    #return r'$f_{{\gamma, {}}}$'.format(fgamma[6:])
def getZetaString(zeta):
    zstr = zeta[4:]
    if zstr[0] == 'D':
        zstr = 'Comp'
    if zstr[-2:] == 'IC':
        return r'$\tilde{{\zeta}}_{{\mathrm{{{}}}}}$'.format(zstr[:-2])
    else:
        return r'$\zeta_{{\mathrm{{{}}}}}$'.format(zstr)

def get_axis_label_from_mod(mod):
    #if 'IC' in mod:
    #    return r'"${}_{\mathrm{var}}$'
    if mod == 'fiducial':
        return mod
    if mod == 'qCritClaeys':
        return r'$q_{{c,\mathrm{Claeys14}}}$' + "\n" + r'($\beta=1.0$)'
    if mod == 'qCritGe20':
        return r'$q_{c,\mathrm{Ge20}}$' + "\n" + r'($\beta=1.0$)'
    if mod == 'qCritGe20IC':
        return r'$\tilde{q}_{c,\mathrm{Ge20}}$' + "\n" + r'($\beta=1.0$)'
    
    # else, must be a zeta
    try:
        zeta, oth = mod.split('_')
        if oth[0] == 'T':
            second_part = "\n" + getTstring(oth)
        else:
            second_part = r'; ' + getBetaString(oth)
    except:
        zeta, beta, fgamma = mod.split('_')
        second_part = r'; ' + getBetaString(beta) + "\n" + getFGammaString(fgamma)
    return getZetaString(zeta) + second_part

#fig, ax = plt.subplots()
#ax.set_xlabel(get_axis_label_from_mod(mods[6]))


# +
def plotMostMods(ax, moddir, fs_axis_label=50, title=None, add_title=False, fs_title=60, addYLabel=True, addXLabel=True, **kwargs):
    
    mods =    [\
        'zetaDefault_betaComp_fgamma0.00',\
        'zetaDefault_betaComp_fgamma1.00',\
        'zetaDefault_beta0.0_fgamma0.00',\
        'zetaDefault_beta0.0_fgamma1.00', \
        'zetaDefault_beta1.0',\
        'qCritClaeys', 'qCritGe20', 'qCritGe20IC']
               #'zetaGe20_beta1.0', 'zetaGe20IC_beta1.0' ,  'zetaGe20_beta0.0_fgamma0.00', 'zetaGe20IC_beta0.0_fgamma0.00', \               'zetaGe20_betaComp_fgamma0.00', 'zetaGe20IC_betaComp_fgamma0.00', 'zetaGe20_betaComp_fgamma1.00', 'zetaGe20IC_betaComp_fgamma1.00', 'zetaGe20_beta0.0_fgamma1.00', 'zetaGe20IC_beta0.0_fgamma1.00'] 
    lw = 6
    
    frac_SMT_end = []
    frac_CEE_end = []
    frac_SMT_mt1 = []
    frac_CEE_mt1 = []
    for ii, modi in enumerate(mods):        
        
        print()
    
        mod = '/{}/{}/'.format(moddir, modi)
        
        ### Extract the post-processed data
        params, stype_masks_colors, state_masks_markers_labels_colors = extractEndStates(mod, printModName=True)#, **kwargs) 
        [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors            

        # Choose whether to plot the mt1_state or the end_state
        for imask, state_masks in enumerate([end_state_masks, mt1_state_masks]):

            data = []
            for stateId in np.unique(state_masks):
                if (stateId == 0):
                    continue # don't include singles
                data.append(np.sum([state_masks == stateId]))        

            fSMT, fCEE = np.cumsum((data/np.sum(data))[1:][::-1]) # Tweak output so that first value is frac of just SMT, second is SMT + CEE surv
            print(fSMT)
            if imask == 0:                
                frac_SMT_end.append(fSMT)
                frac_CEE_end.append(fCEE)
            else:
                frac_SMT_mt1.append(fSMT)
                frac_CEE_mt1.append(fCEE)
                
    
    """
    xPts = np.arange(len(mods), dtype=float) 
    xLastNonGe20 = 7
    xPts[xPts>xLastNonGe20] = xPts[xPts>xLastNonGe20]/2 + xLastNonGe20/2
    """
    xPts = np.arange(len(mods))
    #ax.plot(xPts, frac_CEE_mt1, 'b--', label=r'1st MT: $f_{\mathrm{SMT}} + f_{\mathrm{CEE}}$', lw=lw, alpha=0.5)
    #ax.plot(xPts, frac_SMT_mt1, 'b-',  label=r'1st MT: $f_{\mathrm{SMT}}$', lw=lw)
    #ax.plot(xPts, frac_CEE_end, 'r--', label=r'End state: $f_{\mathrm{SMT}} + f_{\mathrm{CEE}}$', lw=lw, alpha=0.5)
    #ax.plot(xPts, frac_SMT_end, 'r-',  label=r'End state: $f_{\mathrm{SMT}}$', lw=lw)
    ax.plot(xPts, frac_SMT_mt1, '-',  color=colorBlindFriendyPalette[7], label=r'1st MT', lw=lw)
    ax.plot(xPts, frac_SMT_end, '--', color=colorBlindFriendyPalette[2], label=r'End state', lw=lw)
    
    
    ### Setup ticks and labels
    #xMajor = [xPt for xPt in xPts if xPt.is_integer()]
    #xMinor = [xPt for xPt in xPts if not xPt.is_integer()]
    ax.xaxis.set_ticks(xPts)
    #ax.xaxis.set_ticks(xMinor, minor=True )
    

    if addXLabel:
        rotation=50
        fs_modelnames = 34
        xlabels = [get_axis_label_from_mod(mod) for mod in mods] # if not 'IC' in mod] 
        xlabels[0] = xlabels[0] + '*'
        ax.xaxis.set_ticklabels(xlabels, rotation=rotation, fontsize=fs_modelnames)
        #[ax.text(x=xmin, y=-.2, s=r'$\sim$  \;', rotation=rotation, fontsize=fs_modelnames) for xmin in xMinor] #\nwarrow 
        [ax.xaxis.get_majorticklabels()[ii].set_y(-0.1) for ii in range(len(mods))]

    else:
        ax.xaxis.set_ticklabels([])
    ax.tick_params(axis='y', which='major', labelsize=18)    
    ax.tick_params(axis='both', length=10)
    ax.tick_params(axis='both', which='minor', length=10)

    
    
    ###
    ax.set_ylim(-.05, 1.05)
    if addXLabel:
        ax.set_xlabel(r' ', fontsize=fs_axis_label)
    if addYLabel:
        #ax.set_ylabel('Fraction of interacting systems\n undergoing SMT or CEE', fontsize=fs_axis_label)
        ax.set_ylabel(r'$f_{\mathrm{SMT}}$', fontsize=fs_axis_label) #Fraction of systems \n undergoing SMT or CEE

        
    # Setup legend 
    fs_legend = 24
    handles, labels = ax.get_legend_handles_labels()
    #handles.append(mpl.lines.Line2D([0, 0], [0, 0], marker=r'$\sim$', label=r'Ge+20 isentropic variant', color='k', alpha=1.0, markersize=15, lw=0))
    ax.legend(handles=handles, fontsize=fs_legend, framealpha=1, loc='upper left', ncol=1)
    
    
    if title is not None:
        add_title = True
    if add_title:
        if title is None: 
            modcat, moddir, modname = mod.split('/')
            title = "{} :: {} :: {}".format(modcat, moddir, modname)
        ax.set_title(title, fontsize=fs_title)
        
        
    groupbox_color, groupbox_fontsize, groupbox_yval = 'wheat', 20, -0.06
    props = dict(boxstyle='round', facecolor=groupbox_color, alpha=1)
    ax.text(0.03, groupbox_yval, 'x'*43, color=groupbox_color, fontsize=groupbox_fontsize, transform=ax.transAxes, verticalalignment='top', bbox=props)
    ax.text(0.54,  groupbox_yval, 'x'*43, color=groupbox_color, fontsize=groupbox_fontsize, transform=ax.transAxes, verticalalignment='top', bbox=props)
    ax.text(0.18, groupbox_yval, 'Non-conservative', color='k', fontsize=groupbox_fontsize, transform=ax.transAxes, verticalalignment='top')    
    ax.text(0.72, groupbox_yval, 'Conservative', color='k', fontsize=groupbox_fontsize, transform=ax.transAxes, verticalalignment='top')    
    
    
#testing
#fig, ax = plt.subplots(figsize=(20, 8))
#plotMostMods(ax, moddir='uniformIC', title='Traditional Sampler')

# Plot both IC vars side by side
#fig, axes = plt.subplots(ncols=2, figsize=(30, 10))
#plotMostMods(axes[0], moddir='uniformIC', fs_axis_label=36, title='Traditional Sampler', addXLabel=False)
#plotMostMods(axes[1], moddir='moedistef', fs_axis_label=36, title='Moe \& Di Stefano (2017)')

for moddir in ['uniformIC']:#, 'moedistef']:
    fig, ax = plt.subplots(ncols=1, figsize=(15, 10))
    plotMostMods(ax, moddir=moddir, fs_axis_label=36, addXLabel=True)#, title='Traditional Sampler', add_title=False)#, addXLabel=False)
    #plotMostMods(axes[1], moddir='moedistef', fs_axis_label=36, title='Moe \& Di Stefano (2017)')


    fig.tight_layout()
    plt.savefig('../plots/final_plots/mostModsCompared_{}.pdf'.format(moddir), format='pdf', bbox_inches='tight')


# +
def plotMostMods(ax, moddirs, fs_axis_label=50, useEndState=True,  title=None, add_title=False, fs_title=60, addYLabel=True, addXLabel=True, **kwargs):
    
    mods =    [\
        'zetaDefault_betaComp_fgamma0.00',\
        'zetaDefault_betaComp_fgamma1.00',\
        'zetaDefault_beta0.0_fgamma0.00',\
        'zetaDefault_beta0.0_fgamma1.00', \
        'zetaDefault_beta1.0',\
        'qCritClaeys', 'qCritGe20', 'qCritGe20IC']
               #'zetaGe20_beta1.0', 'zetaGe20IC_beta1.0' ,  'zetaGe20_beta0.0_fgamma0.00', 'zetaGe20IC_beta0.0_fgamma0.00', \               'zetaGe20_betaComp_fgamma0.00', 'zetaGe20IC_betaComp_fgamma0.00', 'zetaGe20_betaComp_fgamma1.00', 'zetaGe20IC_betaComp_fgamma1.00', 'zetaGe20_beta0.0_fgamma1.00', 'zetaGe20IC_beta0.0_fgamma1.00'] 
    lw = 6
    
    frac_SMT_end = []
    frac_CEE_end = []
    frac_SMT_mt1 = []
    frac_CEE_mt1 = []
    for ii, modi in enumerate(mods):        
    
        for jj, moddir in enumerate(moddirs):
            
            mod = '/{}/{}/'.format(moddir, modi)
            
            ### Extract the post-processed data
            params, stype_masks_colors, state_masks_markers_labels_colors = extractEndStates(mod, printModName=False, **kwargs) 
            [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors            

            # Choose whether to plot the mt1_state or the end_state
            for imask, state_masks in enumerate([end_state_masks, mt1_state_masks]):

                data = []
                for stateId in np.unique(state_masks):
                    if (stateId == 0):
                        continue # don't include singles
                    data.append(np.sum([state_masks == stateId]))        

                fSMT, fCEE = np.cumsum((data/np.sum(data))[1:][::-1]) # Tweak output so that first value is frac of just SMT, second is SMT + CEE surv
                if imask == 0:                
                    frac_SMT_end.append(fSMT)
                    frac_CEE_end.append(fCEE)
                else:
                    frac_SMT_mt1.append(fSMT)
                    frac_CEE_mt1.append(fCEE)
                
    
    """
    xPts = np.arange(len(mods), dtype=float) 
    xLastNonGe20 = 7
    xPts[xPts>xLastNonGe20] = xPts[xPts>xLastNonGe20]/2 + xLastNonGe20/2
    """
    xPts = np.arange(2*len(mods))
    #ax.plot(xPts, frac_CEE_mt1, 'b--', label=r'1st MT: $f_{\mathrm{SMT}} + f_{\mathrm{CEE}}$', lw=lw, alpha=0.5)
    #ax.plot(xPts, frac_SMT_mt1, 'b-',  label=r'1st MT: $f_{\mathrm{SMT}}$', lw=lw)
    #ax.plot(xPts, frac_CEE_end, 'r--', label=r'End state: $f_{\mathrm{SMT}} + f_{\mathrm{CEE}}$', lw=lw, alpha=0.5)
    #ax.plot(xPts, frac_SMT_end, 'r-',  label=r'End state: $f_{\mathrm{SMT}}$', lw=lw)
    ax.plot(xPts, frac_SMT_mt1, '-',  color=colorBlindFriendyPalette[7], label=r'1st MT', lw=lw)
    ax.plot(xPts, frac_SMT_end, '--', color=colorBlindFriendyPalette[2], label=r'End state', lw=lw)
    
    ### Setup ticks and labels
    #xMajor = [xPt for xPt in xPts if xPt.is_integer()]
    #xMinor = [xPt for xPt in xPts if not xPt.is_integer()]
    ax.xaxis.set_ticks(xPts)
    #ax.xaxis.set_ticks(xMinor, minor=True )
    

    if addXLabel:
        rotation=50
        fs_modelnames = 34
        xlabels = [[get_axis_label_from_mod(mod), r'${}_{\mathrm{[MdS]}}$'][jj] for mod in mods for jj in range(2) ] # if not 'IC' in mod] 
        ax.xaxis.set_ticklabels(xlabels, rotation=rotation, fontsize=fs_modelnames)
        #[ax.text(x=xmin, y=-.2, s=r'$\sim$  \;', rotation=rotation, fontsize=fs_modelnames) for xmin in xMinor] #\nwarrow 

    else:
        ax.xaxis.set_ticklabels([])
    ax.tick_params(axis='y', which='major', labelsize=18)    
    ax.tick_params(axis='both', length=10)
    ax.tick_params(axis='both', which='minor', length=10)

    
    
    ###
    ax.set_ylim(-.05, 1.05)
    if addXLabel:
        ax.set_xlabel(r'Model', fontsize=fs_axis_label)
    if addYLabel:
        #ax.set_ylabel('Fraction of interacting systems\n undergoing SMT or CEE', fontsize=fs_axis_label)
        ax.set_ylabel(r'$f_{\mathrm{SMT}}$', fontsize=fs_axis_label) #Fraction of systems \n undergoing SMT or CEE

        
    # Setup legend 
    fs_legend = 30
    handles, labels = ax.get_legend_handles_labels()
    #handles.append(mpl.lines.Line2D([0, 0], [0, 0], marker=r'$\sim$', label=r'Ge+20 isentropic variant', color='k', alpha=1.0, markersize=15, lw=0))
    ax.legend(handles=handles, fontsize=fs_legend, framealpha=1, loc='upper left', ncol=1)
    
    
    if title is not None:
        add_title = True
    if add_title:
        if title is None: 
            modcat, moddir, modname = mod.split('/')
            title = "{} :: {} :: {}".format(modcat, moddir, modname)
        ax.set_title(title, fontsize=fs_title)
        
#testing
#fig, ax = plt.subplots(figsize=(20, 8))
#plotMostMods(ax, moddir='uniformIC', title='Traditional Sampler')

# Plot both IC vars side by side
fig, ax = plt.subplots(nrows=1, figsize=(25, 10))
plotMostMods(ax, moddirs=['uniformIC','moedistef'], fs_axis_label=36)#, title='Traditional Sampler', addXLabel=False)
#plotMostMods(axes[1], moddir=, fs_axis_label=36, title='Moe \& Di Stefano (2017)')

fig.tight_layout()
plt.savefig('../plots/final_plots/mostModsCompared.pdf', format='pdf', bbox_inches='tight')
# -







if False:
    ### Print latex for appendix plot text generator

    moddir = moddirs[0] # just do one at a time
    mods_list_ordered = [\
     'fiducial',
     'beta0.0',
     'beta0.5',
     'beta1.0',
     'fgamma0.0',
     'fgamma0.25',
     'fgamma0.5',
     'fgamma1.0',
     'fgamma1.0_beta0.0',
     'TbasedConvEnv',
     'zetaGiantInf',
     'zetaMsInf'
    ]


    def full_string(fname):
         return r"""
    \begin{figure}
    \centering
    \includegraphics[width=\columnwidth]{figs/plots_for_appendix/""" \
    + fname + \
    """}
    %\caption{}
    %\label{fig:cdfs}
    \end{figure}
    """


    for mod in mods_list_ordered:
        print(full_string(moddir+'_'+mod+'.png'))









