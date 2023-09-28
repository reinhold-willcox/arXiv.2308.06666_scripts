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

# # Get positions of the lagrange points in terms of the binary separation a, for arbitrary q
#

import numpy as np

# +
s0_s1 = [
    [-1, -1],
    [-1, 1],
    [1, 1],
]

# Ma at 0, Md at -1

def get_alpha(q):
    # q = Ma/Md = M2/M1 
    return 1/(1 + 1/q)

def get_beta(q):
    # q = Ma/Md = M2/M1 
    return 1/(1 + q)
    
    
def get_Lagrange_123(q):
    # q = Ma/Md = M2/M1
    alpha = get_alpha(q)
    
    allroots = []
    for (s0, s1) in s0_s1:
        poly = [
            1,
            3 - alpha,
            3 - 2*alpha,
            1 - s1 - alpha*(1 + s0 - s1),
            -2*alpha*s0,
            -alpha*s0
        ]
        roots = np.roots(poly)
        realroot = roots[np.isreal(roots)].real
        allroots.append(realroot)
    
    L1 = allroots[1]
    if q > 1:
        L2 = allroots[0]
        L3 = allroots[2]
    else:
        L2 = allroots[2]
        L3 = allroots[0]
    return np.array([L1, L2, L3]).flatten()

def get_positions_rel_to_CoM(q):
    L1, L2, L3 = get_Lagrange_123(q)
    posD = -1
    posA = 0
    xCom = -get_beta(q)
    l1l2l3pDpA_relToCom = np.array([L1, L2, L3, posD, posA]) - xCom
    return l1l2l3pDpA_relToCom

#get_Lagrange_123(2)
#get_positions_rel_to_CoM(30)


# +
# Show plot

if False:

    Q = np.logspace(-3, 3, 61)

    cols = ['blue', 'green', 'red']
    labs = ['L1',  'L2', 'L3']

    fig, ax = plt.subplots(figsize=(15, 6))

    for ii, q in enumerate(Q):
        roots = get_Lagrange_123(q)
        for jj, root in enumerate(roots):

            lab = labs[jj]
            col = cols[jj]

            ax.plot(root, q, 'o', color=col, label=lab)

        # CoM
        xCom = -get_beta(q)
        ax.plot(xCom, q, 'y*', label='CoM')

        if ii % 5 == 0:
            ax.plot(-1, q, 'ko', markersize=5*np.sqrt(np.sqrt(1/q))) # M1
            ax.plot(0, q, 'ko', markersize=5*np.sqrt(np.sqrt(q))) # M1

        if ii == 0:
            lgnd = ax.legend(fontsize=16)


    ax.text(x=-1.2, y=1.5, s="M1", fontsize=20)
    ax.text(x=0.1, y=1.5, s="M2", fontsize=20)

    ax.set_yscale('log')
    ax.set_ylabel('q = M2/M1')
    ax.set_xlabel('Position, relative to M2 [a]')
    ax.set_title('Relative positions of CoM, L1, L2, and L3 points for different mass ratios')
    print()
# -

# # Conclusions:
#
# ### For q != 1, the CoM is closer to the larger component, but the L1 point is closer to the smaller one
#
# ### The CoM approaches its nearest component faster than the L points 

# +

# Want to get separation from CoM to L2 point, see if its consistently close to 1.2a

if False:
    Q = np.logspace(-3, 3, 101)
    #print(Q)

    fig, axes = plt.subplots(ncols=2, figsize=(15, 8))

    Sep = []

    for q in Q:
        sep = get_sep_L2_CoM(q)
        Sep.append(sep)
        
    for ii, ax in enumerate(axes):

        ax.plot(Q, Sep)
        ax.set_xscale('log')
        ax.set_xlabel('q = M2/M1')

        if ii == 0:
            ax.set_ylabel(r'$d_{\mathrm{L2}}$ [a]')
            ax.set_title('Distance of L2 point from CoM')
        else:
            ax.set_xlim((.1, 10))
            ax.set_title('Same plot, zoomed in')
# -




