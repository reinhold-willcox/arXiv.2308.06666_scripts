# ---
# jupyter:
#   jupytext:
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

import cmasher as cmr
from Lagrange_points import get_positions_rel_to_CoM

# +
### For paper, plot beta vs q for a MS -> MS MT, for a variety of C's

fig, ax = plt.subplots(figsize=(10, 6))

q = np.linspace(0, 1, 10001)
q3_5 = np.power(q, 3.5)
print(q3_5)
C = [1, 3, 10]
ls = ['r-', 'g--', 'b:']

for ii in range(len(C)):
    beta = np.minimum(C[ii] * q3_5, 1)
    ax.plot(q, beta, ls[ii], label='C='+str(C[ii]), lw=3)
    
ax.legend(fontsize=16)
fs=30
ax.set_xlabel('q', fontsize=fs)
ax.set_ylabel(r'$\beta$', fontsize=fs)

plt.savefig('../plots/final_plots/beta_vs_q.pdf', format='pdf', bbox_inches='tight')
# -







# +
q = np.linspace(0, 1, 101)

def dlnrl_dlnq(q):
    q13 = np.power(q, 1/3)
    q23 = q13*q13
    num1 = (0.17/q13)*(0.6*q23 + np.log(1+q13))
    num2 = (0.49*q23)*(0.4/q13 + 1/(3*q23*(1+q13)))
    den = (0.6*q23 + np.log(1+q13))
    return (num1-num2)/(den*den)



# -

fig, ax = plt.subplots(figsize=(10, 4))
fontsize = 30
ax.plot(q, dlnrl_dlnq(q))
ax.set_ylabel(r'$\frac{d\ln(r_L)}{d\ln(q)}$', fontsize=fontsize)
ax.set_xlabel(r'$q$', fontsize=fontsize)
ax.set_xscale('log')
plt.savefig('../plots/dlnrLdlnq_vs_q.eps', format='eps', bbox_inches='tight')






# +
### Plot Soberman zeta_SPH as a function of core mass

def zeta_SPH(m):
    # m is m_core/m_total
    
    zeta_sph1 = (2/3)*(m/(1-m)) - (1/3)*((1-m)/(1+2*m))
    
    one_minus_m_neg6 = np.power(1-m, -6)    
    rhs = -0.03*m + 0.2* (m/(1+one_minus_m_neg6))

    return zeta_sph1 + rhs

m_frac = np.linspace(0, 1, 101)[:-1]
print(m_frac)

# +
fig, ax = plt.subplots(figsize=(15, 8))

ax.plot(m_frac, zeta_SPH(m_frac))
ax.set_ylim(-2, 10)
fs = 30
ax.set_ylabel(r'$\zeta_{{SPH}}$', fontsize=fs)
ax.set_xlabel(r'Core mass fraction $M_c/M$', fontsize=fs)

plt.savefig('../plots/zeta_sph.eps', format='eps', bbox_inches='tight')
# -





# +
# Plot zeta_a for a grid of q, beta, and a_gamma/a (where a_gamma is the point of detachment)


def get_zeta_a(beta, q, a_rat):
    # uses Ma/Md
    first_term =  1 - beta/q
    second_term_pt1 = (1-beta)*(1/(1+q))
    second_term_pt2 = a_rat*a_rat * (1+q)*(1+q)/q + 1/2
    second_term = second_term_pt1 * second_term_pt2.reshape((-1, 1))
    return -2 *(first_term - second_term)

def get_drL_dq(q):
    # uses Md/Ma
    q13 = np.cbrt(q)
    first_term = q13/3
    second_pt1 = 2/q13
    second_pt2_num = 1.2*q13 + 1/(1+q13)
    second_pt2_den = 0.6*q13*q13 + np.log(1+q13)
    return first_term *(second_pt1 - second_pt2_num/second_pt2_den)

def get_zeta_q(beta, q):
    # uses Md/Ma
    return 1 - q*beta


def get_zeta_L(beta, q, a_rat):
    # q is Ma/Md
    zeta_a = get_zeta_a(beta, q, a_rat)
    #print(zeta_a)
    drL_on_q = get_drL_dq(1/q)
    zeta_q = get_zeta_q(beta, 1/q)
    return zeta_a + drL_on_q * zeta_q
# +
# try again, this time with my own definitions


def calc_zeta_a(beta, q, a_rat):
    # uses Ma/Md
    first_term =  1 - beta/q
    second_term_pt1 = (1-beta)*(1/(1+q))
    second_term_pt2 = a_rat*a_rat * (1+q)*(1+q)/q + 1/2
    second_term = second_term_pt1 * second_term_pt2.reshape((-1, 1))
    return -2 *(first_term - second_term)


def calc_drL_dq(q):
    # uses Ma/Md
    q13 = np.cbrt(q)
    prefactor = -q13*q13/3
    first = 2/(q13*q13)
    top = 1.2/(q13*q13) + 1/(1+q13)
    bottom = 0.6 + q13*q13*np.log(1+1/q13)
    return prefactor* (first - top/bottom)

def calc_drL_dq2(q):
    # uses Ma/Md
    q13 = np.cbrt(q)
    qn13 = 1/q13
    top = 1.2 + 1/(qn13 + qn13*qn13)
    bottom = 0.6 + q13*q13*np.log(1+qn13)
    return -1/3 * (2 - top/bottom)
    
    
def calc_zeta_q(beta, q):
    # uses Ma/Md
    return beta/q - 1


# todo
def calc_zeta_L(beta, q, a_rat, test=False):
    # q is Ma/Md
    zeta_a = calc_zeta_a(beta, q, a_rat)
    drL_on_q = calc_drL_dq(q)
    if test:
        drL_on_q = calc_drL_dq2(q)
    zeta_q = calc_zeta_q(beta, q)
    return zeta_a + drL_on_q * zeta_q


# -

#testing
if True:
    beta = 0.5
    a_rat = .8

    for q in np.logspace(-1, 1, 21):
        #print(q)
        z1 = get_zeta_L(beta, q, a_rat)
        z2 = calc_zeta_L(beta, q, a_rat)
        z3 = calc_zeta_L(beta, q, a_rat, test=True)
        print(z1)
        print(z2)
        print(z3)
        print()








# +
# Plot roche lobe response zeta_RL to mass loss

Q = np.linspace(0, 10, 11)[1:]/10 
beta = np.linspace(0, 1, 201)
a_rat = np.linspace(0, 1.3, 281) 

# L2 position calculated using lagrange solver in other notebook, this is for q = np.linspace(1, 1, 10)
L2_position = [1.2560829084936194, 1.2714101144108156, 1.2684022053216617, 1.259666516797683, 1.2490473888803288, 1.2380383666032049, 1.227274369630634, 1.2170259499466047, 1.2073958135032674, 1.19840614455492]


labY = [.1, .6]
labX = [.6, .7, .8, .9, 1.0]

fs_nums = 18
fs_labels = 30

vmin, vmax = -2, 12
lower = 0.5*(1+vmin/vmax)
colormap = 'seismic' #'cmr.fusion' #'nipy_spectral' # 'cmr.pride' #'seismic' #'hsv'
#cmap = plt.get_cmap('seismic')(np.linspace(lower, 1.0, 100))
cmap = cmr.get_sub_cmap(colormap, lower, 1.0, N=3*(vmax-vmin))


fig, axes = plt.subplots(ncols=5, nrows=2, figsize=(25, 10))
axs = axes.flatten()

for q, ax, l2 in zip(Q, axs, L2_position):
    zeta_L = get_zeta_L(beta, q, a_rat[::-1])
    extent = (beta[0], beta[-1], a_rat[0], a_rat[-1])
    #cmap = 'seismic'
    #vmin, vmax = -6, 6
    im = ax.imshow(zeta_L, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent) #(left, right, bottom, top), optional
    ax.set_title(r"$q = {{}}${:.1f}".format(q), fontsize=fs_labels)
    if q in labY:
        ax.set_ylabel(r'$a_\gamma / a$', fontsize=fs_labels)
    if q in labX:
        ax.set_xlabel(r'$\beta$', fontsize=fs_labels)
        
    ax.tick_params(axis='both', labelsize=fs_nums)

    ax.set_ylim(0.0, l2)
    a_acc_over_a = 1/(1+q)
    ax.axhline(y=a_acc_over_a, color='fuchsia', lw=3)
    
    print(np.min(zeta_L), np.max(zeta_L))
    
fig.subplots_adjust(right=0.85)
#cbar_ax = fig.add_axes([1.03, 0.15, 0.04, 0.7])
cbar_ax = fig.add_axes([0.1, -0.15, .8, 0.1])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label(r'$\zeta_{L} = d\ln(R_{RL})/d\ln(M_d)$', fontsize=fs_labels)

fig.tight_layout()

#plt.savefig('../plots/final_plots/zetaRL_heatmap.pdf', format='pdf', bbox_inches='tight')


# +
# Plot roche lobe response zeta_RL to mass loss

#Q = np.linspace(0, 10, 11)[1:]/10 
Q1 = np.linspace(0.2, 1, 5, endpoint=True)
Q2 = np.array([1.5, 2, 3, 4, 5])
Q = np.append(Q1, Q2)

beta = np.linspace(0, 1, 201)
a_rat = np.linspace(0, 1.3, 281) 


labY = [0, 5]
labX = [5, 6, 7, 8, 9]

fs_nums = 18
fs_labels = 30

vmin, vmax = -2, 10
lower = 0.5*(1+vmin/vmax)
colormap = 'seismic' #'cmr.fusion' #'nipy_spectral' # 'cmr.pride' #'seismic' #'hsv'
#cmap = plt.get_cmap('seismic')(np.linspace(lower, 1.0, 100))
cmap = cmr.get_sub_cmap(colormap, lower, 1.0, N=3*(vmax-vmin))


fig, axes = plt.subplots(ncols=5, nrows=2, figsize=(25, 10))
axs = axes.flatten()

for ii, ax in enumerate(axs):
    q = Q[ii]
    zeta_L = get_zeta_L(beta, q, a_rat[::-1])
    extent = (beta[0], beta[-1], a_rat[0], a_rat[-1])
    #cmap = 'seismic'
    #vmin, vmax = -6, 6
    im = ax.imshow(zeta_L, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent) #(left, right, bottom, top), optional
    ax.set_title(r"$q = {{}}${:.1f}".format(q), fontsize=fs_labels, y=1.02)
    if ii in labY:
        ax.set_ylabel(r'$a_\gamma / a$', fontsize=fs_labels)
    if ii in labX:
        ax.set_xlabel(r'$\beta$', fontsize=fs_labels)
        
    ax.tick_params(axis='both', labelsize=fs_nums)

    l1l2l3pDpA_relToCom = get_positions_rel_to_CoM(q)
    _, l2, l3, _, _ = l1l2l3pDpA_relToCom
    ax.set_ylim(0.0, max(l2, l3))
    a_acc_over_a = 1/(1+q)
    ax.axhline(y=a_acc_over_a, color='fuchsia', lw=3)
    
    print(np.min(zeta_L), np.max(zeta_L))
    
fig.subplots_adjust(right=0.85)
#cbar_ax = fig.add_axes([1.03, 0.15, 0.04, 0.7])
cbar_ax = fig.add_axes([0.1, -0.15, .8, 0.1])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label(r'$\zeta_{L} = d\ln(R_{RL})/d\ln(M_d)$', fontsize=fs_labels)

fig.tight_layout()

plt.savefig('../plots/final_plots/zetaRL_heatmap.pdf', format='pdf', bbox_inches='tight')


# +
# Plot roche lobe response zeta_RL to mass loss
#Q = [30]

Q1 = np.linspace(0.2, 1, 5, endpoint=True)
Q2 = np.array([1.5, 2, 3, 4, 5])
Q = np.append(Q1, Q2)
Q = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 3.0, 5.0]
#Q = [0.03, 0.1, 0.3, 0.6, 0.8, 1.0, 1.5, 5.0, 10, 30]
beta = np.linspace(0, 1, 401)
a_rat = np.linspace(-1.3, 1.3, 271) 

labY = [0, 5]
labX = [5, 6, 7, 8, 9]

fs_nums = 18
fs_labels = 30

vmin, vmax =-2, 12
lower = 0.5*(1+vmin/vmax)
colormap = 'seismic' #'cmr.fusion' #'nipy_spectral' # 'cmr.pride' #'seismic' #'hsv'
#cmap = plt.get_cmap('seismic')(np.linspace(lower, 1.0, 100))
cmap = cmr.get_sub_cmap(colormap, lower, 1.0, N=(vmax-vmin)*4)



fig, axes = plt.subplots(ncols=5, nrows=2, figsize=(25, 13))
axs = axes.flatten()

for ii, ax in enumerate(axs):
    q = Q[ii]
    l1l2l3pDpA_relToCom = get_positions_rel_to_CoM(q)
    l1, l2, l3, pD, pA = l1l2l3pDpA_relToCom
    if l3 > l2:
        l2, l3 = l3, l2

    zeta_L = get_zeta_L(beta, q, a_rat[::-1])
    #extent = (beta[0], beta[-1], a_rat[0], a_rat[-1])#(left, right, bottom, top), optional
    extent = (beta[0], beta[-1], l3, l2)#(left, right, bottom, top), optional
    im = ax.imshow(zeta_L, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent, aspect=0.75) 

    # Set title and axes labels
    ax.set_title(r"$q = {{}}${:.1f}".format(q), fontsize=fs_labels, y=1.02)
    if ii in labY:
        ax.set_ylabel(r'$a_\gamma / a$', fontsize=fs_labels)
    if ii in labX:
        ax.set_xlabel(r'$\beta$', fontsize=fs_labels)

    # Set ticks
    ax.tick_params(axis='both', labelsize=fs_nums)
    ax.set_xticks([0, 0.5, 1.0])
    #ax.set_yticks()

    # Set the text on the right
    xtext=1.05
    
    if l2 > l3:
        yNear = l3
        yFar  = l2
        ax.text(x=xtext, y=yNear, s='L3', fontsize=fs_nums)
        ax.text(x=xtext, y=yFar-.1,  s='L2', fontsize=fs_nums)
    else:
        yNear = l2
        yFar  = l3
        ax.text(x=xtext, y=yNear, s='L2', fontsize=fs_nums)
        ax.text(x=xtext, y=yFar-.1,  s='L3', fontsize=fs_nums)        
    
    for jj, pp in enumerate([l1, pD, pA]):
        col = ['fuchsia', 'y','g'][jj]
        lab = ['L1', r'$M_D$', r'$M_A$'][jj]
        ax.axhline(y=pp, color=col, lw=3)
        ax.text(x=xtext, y=pp, s=lab, fontsize=fs_nums)
    
    print(np.min(zeta_L), np.max(zeta_L))
    
cbar_ax = fig.add_axes([0.1, -0.15, .8, 0.1])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label(r'$\zeta_{L} = d\ln(R_{RL})/d\ln(M_d)$', fontsize=fs_labels)
cbar.set_ticks(list(np.linspace(-2, 12, 15)))#, labelsize=fs_nums)
cbar.ax.tick_params(labelsize=fs_nums)

fig.tight_layout()
plt.savefig('../plots/final_plots/zetaRL_heatmap.pdf', format='pdf', bbox_inches='tight')

# +
# Adelle system

# Plot roche lobe response zeta_RL to mass loss
Q = [30]
beta = np.linspace(0, 1, 401)
a_rat = np.linspace(-1.3, 1.3, 271) 

labY = [0, 5]
labX = [5, 6, 7, 8, 9]

fs_nums = 18
fs_labels = 30

vmin, vmax =-2, 12
lower = 0.5*(1+vmin/vmax)
colormap = 'seismic' #'cmr.fusion' #'nipy_spectral' # 'cmr.pride' #'seismic' #'hsv'
#cmap = plt.get_cmap('seismic')(np.linspace(lower, 1.0, 100))
cmap = cmr.get_sub_cmap(colormap, lower, 1.0, N=(vmax-vmin)*4)



fig, ax = plt.subplots(figsize=(12, 15))

q = Q[0]
l1l2l3pDpA_relToCom = get_positions_rel_to_CoM(q)
l1, l2, l3, pD, pA = l1l2l3pDpA_relToCom

zeta_L = get_zeta_L(beta, q, a_rat[::-1])
extent = (beta[0], beta[-1], a_rat[0], a_rat[-1])#(left, right, bottom, top), optional
im = ax.imshow(zeta_L, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent, aspect=0.75) 

# Set title and axes labels
ax.set_title(r"$q = {{}}${:.1f}".format(q), fontsize=fs_labels, y=1.02)
if ii in labY:
    ax.set_ylabel(r'$a_\gamma / a$', fontsize=fs_labels)
if ii in labX:
    ax.set_xlabel(r'$\beta$', fontsize=fs_labels)

# Set ticks
ax.tick_params(axis='both', labelsize=fs_nums)
ax.set_xticks([0, 0.5, 1.0])
#ax.set_yticks()

# Set the text on the right
xtext=1.05
if l2 > l3:
    yNear = l3
    yFar  = l2
    ax.text(x=xtext, y=yNear, s='L3', fontsize=fs_nums)
    ax.text(x=xtext, y=yFar-.1,  s='L2', fontsize=fs_nums)
else:
    yNear = l2
    yFar  = l3
    ax.text(x=xtext, y=yNear, s='L2', fontsize=fs_nums)
    ax.text(x=xtext, y=yFar-.1,  s='L3', fontsize=fs_nums)        
ax.set_ylim(yNear, yFar)
for jj, pp in enumerate([l1, pD, pA]):
    col = ['fuchsia', 'y','g'][jj]
    lab = ['L1', r'$M_D$', r'$M_A$'][jj]
    ax.axhline(y=pp, color=col, lw=3)
    ax.text(x=xtext, y=pp, s=lab, fontsize=fs_nums)

print(np.min(zeta_L), np.max(zeta_L))
    
cbar_ax = fig.add_axes([0.1, -0.15, .8, 0.1])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label(r'$\zeta_{L} = d\ln(R_{RL})/d\ln(M_d)$', fontsize=fs_labels)
cbar.set_ticks(list(np.linspace(-2, 12, 15)), fontsize=fs_nums)

fig.tight_layout()

# -





# +
# For a few fixed q values (and holding a fixed), calculate 





def evolve_RL_from_MT(beta, fgamma, q0, dM_don = .01, nTimesteps = 100):
    Md = [1, 1]
    Rl = [1, 1]
    Ma = [q0, q0]
    dM_acc = beta*dM_don
    
    while (Rl[-1] > 0.1) and (Rl[-1] <= Rl[-2]):    
        Md_post = Md[-1] - dM_don
        Ma_post = Ma[-1] + dM_acc
        
        q = Ma[-1]/Md[-1]
        arat_fg0 = 1/(1+q) # this is fgamma=0
        arat_fg1 = 1.2
        arat = np.array([arat_fg0 + fgamma*(arat_fg1-arat_fg0)])

        zeta  = get_zeta_L(beta, q, arat).flatten()[0]
        Rl_post = Rl[-1] - (Md[-1]/Rl[-1])*zeta*dM_don
        
        Md.append(Md_post)
        Ma.append(Ma_post)
        Rl.append(Rl_post)
        
    return Md, Ma, Rl

# testing
Md, Ma, Rl = evolve_RL_from_MT(1, 0, 0.5)
#Rl




beta = [0, 0.5, 1]
fgamma = [0, 0.5, 1]
Q0 = np.linspace(0, 1, 11, endpoint=True)[1:] # the initial q value
print(Q0)

fig, axes = plt.subplots(ncols=3, figsize=(20, 8))

for ii in range(3):
    for jj in range(3):
        for kk in range(len(Q0)):
            ax = axes[ii] # diff fg
            
            bt = beta[ii]
            fg = fgamma[jj]
            q0 = Q0[kk]
            
            col = ['r', 'g', 'k'][jj] # color on beta
            #ls = ['-', '--', ':'][kk]


            Md, Ma, Rl = evolve_RL_from_MT(bt, fg, q0)
            ax.plot(Md, np.log10(Rl), c=col)#, ls=ls)
