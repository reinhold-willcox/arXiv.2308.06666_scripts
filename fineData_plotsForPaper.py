# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # This plot combines and formats the plots made in individual_plotsForPaper
from stable_mt_utils_fix import extractEndStates

# +
colorBlindFriendyPalette = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

fig, ax = plt.subplots()

for ii in range(len(colorBlindFriendyPalette)):
    ax.plot([0, 1], [ii, ii], color=colorBlindFriendyPalette[ii], lw=8)


# -

def plotGammaVsBeta(ax, moddir, fs_axis_label=30, yUpp=None, altPrimaryDir=False, includeFCEE=False, useEndState=True, zeta='Default', addZeta=False, title=None, add_title=False, addExplanation=False, corner_letter=None, addYLabel=True, addXLabel=True, fs_title=40, **kwargs):
    
    betas = ['1.0', '0.5', '0.0', 'Comp']#, 'CompC03', 'CompC30']
    #colors = ['red', 'green', 'blue', 'k', 'orange', 'pink']
    colors = [colorBlindFriendyPalette[7], colorBlindFriendyPalette[5], colorBlindFriendyPalette[2], 'k']
    linestyles = ['--', ':', '-.', '-'] #, ':', '-.', '-']0.00
    if altPrimaryDir:
        fgs = ['0.00', '0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90', '1.00'] #
    else:
        fgs =  [ '0.00', '0.25', '0.50', '0.75', '1.00'] #
    fg_floats = [float(fg) for fg in fgs]
    lw = 4

    grabDataFromH5 = True
    #with open('finely_sampled_data.dat', 'a') as f:
    for beta, col, ls in zip(betas, colors, linestyles):
        
        print()
        
        if grabDataFromH5:
            frac_SMT = []
            frac_CEE = []

            for fg in fgs:
                mod = '{}/zeta{}_beta{}_fgamma{}/'.format(moddir, zeta, beta, fg)
                if beta == '1.0':
                    mod = '{}/zeta{}_beta{}/'.format(moddir, zeta, beta, fg)
                ### Extract the post-processed data
                print(mod)
                params, stype_masks_colors, state_masks_markers_labels_colors = extractEndStates(mod, printModName=True)#, altPrimaryDir=altPrimaryDir, **kwargs) 
                [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors] = state_masks_markers_labels_colors            

                # Choose whether to plot the mt1_state or the end_state
                if useEndState:
                    state_masks = end_state_masks 
                else:
                    state_masks = mt1_state_masks 

                data = []
                for stateId in np.unique(state_masks):
                    if (stateId == 0):
                        continue # don't include singles
                    data.append(np.sum([state_masks == stateId]))

                fSMT, fCEE = np.cumsum((data/np.sum(data))[1:][::-1]) # Tweak output so that first value is frac of just SMT, second is SMT + CEE surv
                print(fSMT)
                frac_SMT.append(fSMT)
                frac_CEE.append(fCEE)
                
        else:
            # read in from .dat file
            fname = '../data/finely_sampled_data_{}.dat'.format('end' if useEndState else 'mt1')
            mod = '{}/zeta{}_beta{}_fgamma{}/'.format(moddir, zeta, beta, '1.00')
            if beta == '1.0':
                mod = '{}/zeta{}_beta{}/'.format(moddir, zeta, beta)

            with open(fname, 'r') as f:
                lines = f.read().splitlines()
            index = lines.index(mod)
            frac_CEE = np.array(lines[index+1][1:-1].split(', '), dtype=float)
            frac_SMT = np.array(lines[index+2][1:-1].split(', '), dtype=float)
  

        ax.plot(fg_floats, frac_SMT, lw=lw, ls=ls, color=col, alpha=1)#, label='SMT : ' + beta)
        if includeFCEE:
            ax.plot(fg_floats, frac_CEE, lw=lw, ls=linestyles[1], color=col, alpha=0.5)#, label='CEE : ' + beta)
        #f.write('\n' + mod + '\n')
        #f.write('[' + ', '.join(["{:.2f}".format(frac) for frac in frac_CEE]) + ']\n')
        #f.write('[' + ', '.join(["{:.2f}".format(frac) for frac in frac_SMT]) + ']\n')
        #print(mod)
        #print(frac_CEE)
        #print(frac_SMT)
        #print()
    
    
        
    
    # Set label names
    if addXLabel:
        ax.set_xlabel(r'Specific angular momentum, $f_\gamma$', fontsize=fs_axis_label)
    if addYLabel:
        if includeFCEE:
            ax.set_ylabel(r'$f_{\mathrm{SMT}}, f_{\mathrm{CEE}}$', fontsize=fs_axis_label) #Fraction of systems \n undergoing SMT or CEE
        else:
            ax.set_ylabel(r'$f_{\mathrm{SMT}}$', fontsize=fs_axis_label) #Fraction of systems \n undergoing SMT or CEE
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.yaxis.set_ticks(np.linspace(0, 1, 11))

    
    # Set y limits
    if yUpp is None:
        if zeta == 'Default':
            ax.set_ylim(-.05, .55)
        else:
            ax.set_ylim(-.05, .7)
    else:
        if yUpp < 0.5:
            #ax.yaxis.set_ticks(np.linspace(0, 1, 21), which='minor')
            minor_locator = mpl.ticker.AutoMinorLocator(2)
            ax.yaxis.set_minor_locator(minor_locator)
            ax.grid(which='minor')
        ax.set_ylim(-.04, yUpp)
    


    
    # Add custom legend 
    def beta_string(beta):
        return r'$\beta_{{Comp}}$' if beta == 'Comp' else r'$\beta=$'+beta
    custom_legend = []
    label_prefix = '' #'Pre-SN: ' if useEndState else '1st MT: '
    if includeFCEE:
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=label_prefix + r'$f_{\mathrm{SMT}} + f_{\mathrm{CEE}}$', color='k', alpha=0.7, lw=2, linestyle=linestyles[1]))
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=label_prefix + r'$f_{\mathrm{SMT}}$', color='k', alpha=1.0, lw=2, linestyle=linestyles[0]))
    for ii in range(len(betas)):#[::-1]:
        custom_legend.append(mpl.lines.Line2D([0, 0], [0, 0], label=beta_string(betas[ii]), color=colors[ii], alpha=1.0, lw=3.5, linestyle=linestyles[ii]))
    ax.legend(handles=custom_legend, fontsize=30, framealpha=1, loc='lower left', ncol=1)
    
    
    # Add ABCD corner letters
    if corner_letter is not None:
        corner_letter = ['a', 'b', 'c', 'd'][corner_letter]
        corner_string = "{})".format(corner_letter)
        ax.text(x=.02, y=ax.get_ylim()[1]*0.93, s=corner_string, fontsize=24)
        #print('yes')
    
    # Set title
    if title is not None:
        add_title = True
    if add_title:
        if title is None: 
            modcat, moddir, modname = mod.split('/')
            title = "{} :: {} :: {}".format(modcat, moddir, modname)
        ax.set_title(title, fontsize=fs_title)
        
    # Add Zetas
    if addZeta:
        zetaStr = r'$\zeta_\mathrm{{ {} }}$'.format(zeta)
        ax.text(x=-0.25, y=sum(ax.get_ylim())*0.4, s=zetaStr, fontsize=fs_title*1.2, rotation=90)
        
    # Add explanation
    if addExplanation:
        fs_explanation = 18
        x_explanation = .02
        x_lines = .7
        ax.annotate(text='', xy=(x_lines,.15), xytext=(x_lines,.285), arrowprops=dict(arrowstyle='<->', color='red'))
        ax.annotate(text='', xy=(x_lines,0), xytext=(x_lines,.15), arrowprops=dict(arrowstyle='<->', color='red'))
        ax.annotate(text='', xy=(x_lines,0.23), arrowprops=dict(arrowstyle='<-'), xytext=(.35,.13), fontsize=12)
        ax.annotate(text='', xy=(x_lines,0.08), arrowprops=dict(arrowstyle='<-'), xytext=(.3,.07), fontsize=12)
        ax.annotate(text=r'$\sim13\%$ CEE survivors',                    xytext=(x_explanation,.12), fontsize=fs_explanation, xy=(.3,0.18))
        ax.annotate(text=r'$\sim15\%$ Only SMT',                         xytext=(x_explanation,.06), fontsize=fs_explanation, xy=(.3,0.05))
        ax.annotate(text=r'Thus $\sim72\%$ merge during the first MT',   xytext=(x_explanation,.01), fontsize=fs_explanation, xy=(.3,0.05))

# +
fig, ax = plt.subplots(figsize=(10, 8))

plotGammaVsBeta(ax, moddir='uniformIC', yUpp=.7, zeta='Default', useEndState=False, addZeta=False, addXLabel=True, addYLabel=True, title=None, addExplanation=False, corner_letter=None)

# -



# +
# Plot files individually

#for ii in range(len(axs)):


#useEndState = True # do both!
moddirs = ['uniformIC', 'moedistef']
titles = ['Traditional Sampler', 'Moe \& Di Stefano (2017)', None, None]
addX = [False, False, True, True]
addY = [True, False, True, False]
zetas = ['Default', 'Default', 'Ge20', 'Ge20']
addZetas = [True, False, True, False]

for useEndState in [True, False]:
    for ii in range(2):

        fig, ax = plt.subplots(figsize=(10, 8))
        
        yUpp = 0.3 if useEndState else 0.6
        addExp = [not useEndState, False, False, False]

        moddir=moddirs[ii%2]
        #zeta=zetas[ii]
        #plotGammaVsBeta(ax, moddir=moddirs[ii%2], yUpp=yUpp, zeta=zetas[ii], useEndState=useEndState, addZeta=addZetas[ii], addXLabel=addX[ii], addYLabel=addY[ii], title=titles[ii], addExplanation=addExp[ii], corner_letter=None)
        plotGammaVsBeta(ax, moddir=moddirs[ii%2], yUpp=yUpp, zeta='Default', useEndState=useEndState, addZeta=False, addXLabel=True, addYLabel=True, title=None, addExplanation=False, corner_letter=None)

        #fig.tight_layout()

        endStateStr = 'end' if useEndState else 'mt1'

        plt.savefig('../plots/final_plots/gammaVsBeta_{}_{}.pdf'.format(endStateStr, moddir), format='pdf', bbox_inches='tight')
# -



