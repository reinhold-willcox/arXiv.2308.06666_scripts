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

import numpy as np
import h5py as h5

def get_h5(mod, altPrimaryDir=False):
    if altPrimaryDir:
        return h5.File('../data/frivolous_extra_gammas/{}/COMPAS_Output.h5'.format(mod), 'r')
    else:
        return h5.File('../data/masterOutput/{}/COMPAS_Output.h5'.format(mod), 'r')

def checkMasksMutuallyExclusive(masks):
    # Runs assert as the check - if this returns with no fuss, the masks are exclusive
    if len(masks) < 2:
        return

    ## verify the state masks
    for mm in range(len(masks)):
        for nn in range(mm):
            msk1 = masks[mm]
            msk2 = masks[nn]
            assert len(msk1) == len(msk2)
            assert ~np.any(msk1 & msk2)



def extractEndStates(mod, printModName=True, altPrimaryDir=False):

    # testing
    #mod = 'uniformIC/fiducial'
    #altPrimaryDir = False
    #printModName = False

    myf = get_h5(mod, altPrimaryDir)
    if printModName:
        print(mod)

    SPs = myf['BSE_System_Parameters']
    MTs = myf['BSE_RLOF']
    CEs = myf['BSE_Common_Envelopes']

    # Get parameters
    spSeeds = SPs['SEED'][()]
    mtSeeds = MTs['SEED'][()]
    ceSeeds = CEs['SEED'][()]
    m1 = SPs['Mass@ZAMS(1)'][()]
    m2 = SPs['Mass@ZAMS(2)'][()]
    mT = m1+m2
    st1_f = SPs['Stellar_Type(1)'][()]
    st2_f = SPs['Stellar_Type(2)'][()]
    aZams = SPs['SemiMajorAxis@ZAMS'][()] # AU
    eZams = SPs['Eccentricity@ZAMS'][()]    

    st1_mt= MTs['Stellar_Type(1)<MT'][()]
    st2_mt= MTs['Stellar_Type(2)<MT'][()]
    mtMerger = MTs['Merger'][()] == 1
    mtCEE    = MTs['CEE>MT'][()] == 1
    isRlof1 = MTs['RLOF(1)>MT'][()] == 1
    isRlof2 = MTs['RLOF(2)>MT'][()] == 1
    isRlof1Prev = MTs['RLOF(1)<MT'][()] == 1
    isRlof2Prev = MTs['RLOF(2)<MT'][()] == 1

    st1_ce = CEs['Stellar_Type(1)<CE'][()]
    st2_ce = CEs['Stellar_Type(2)<CE'][()]

    # Derived quantities
    q = m2/m1
    rP = aZams*(1-eZams)



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
    maskPreSN = ((st1_mt < 13) | (st1_mt == 16)) & ((st2_mt < 13) | (st2_mt == 16))  # 16 for CHE which are MS stars for this purpose
    
    mask1stMtFirstTimestep = mask1stMT & maskPreSN
    mask1stMtLastTimestep  = (mask1stMtFirstTimestep & ~maskMtStartOfContinuation) 
    # need to add into it, the maskMtend ofcontinuation but only for the ones that apply for the 1st MT
    mask1stMtLastTimestep[maskMtEndOfContinuation] = mask1stMtFirstTimestep[maskMtStartOfContinuation]    


    maskStable1stMt       = mask1stMtLastTimestep & ~mtCEE & ~mtMerger 
    mask1stMerger         = mask1stMtLastTimestep & mtMerger           
    mask1stMergerStable   = mask1stMtLastTimestep & mtMerger & ~mtCEE  
    mask1stMergerUnStable = mask1stMtLastTimestep & mtMerger & mtCEE   



    st1mt_1stPreSnMt = st1_mt[mask1stMtFirstTimestep]
    st1mt_1stPreSnMt[st1mt_1stPreSnMt == 16] = 1 # Reset stype 16 (CHE) to 1
    # Need to remap these types back onto a SP-like array
    sdsMt_1stPreSnMt = mtSeeds[mask1stMtFirstTimestep]
    spMask_1stPreSnMt = np.in1d(spSeeds, sdsMt_1stPreSnMt)
    st1sp_1stPreSnMt = np.zeros_like(spSeeds)                # initialize all stypes to 0, to hide these away
    st1sp_1stPreSnMt[spMask_1stPreSnMt] = st1mt_1stPreSnMt   # update the interactors to the correct stype    
    stype_masks = st1sp_1stPreSnMt


    

    ### Extract the end state after the 1st MT (isolated, only SMT, CEE, and merger) 

    maskMtMerger1stMT  = mask1stMtLastTimestep & mtMerger
    maskMtCeeSurv1stMT = mask1stMtLastTimestep & mtCEE & ~maskMtMerger1stMT
    maskMtSMT1stMT     = mask1stMtLastTimestep & ~maskMtMerger1stMT & ~maskMtCeeSurv1stMT
    maskSpMerger1stMT  = np.in1d(spSeeds, mtSeeds[maskMtMerger1stMT]) 
    maskSpCeeSurv1stMT = np.in1d(spSeeds, mtSeeds[maskMtCeeSurv1stMT])
    maskSpSMT1stMT     = np.in1d(spSeeds, mtSeeds[maskMtSMT1stMT])    

    mt1_state_mask_list = [maskSpMerger1stMT , maskSpCeeSurv1stMT, maskSpSMT1stMT]
    checkMasksMutuallyExclusive(mt1_state_mask_list) # check that things look good
    mt1_state_masks = np.zeros_like(spSeeds)
    for ii, mask in enumerate(mt1_state_mask_list):
        maskIdNumber = ii + 1
        mt1_state_masks[mask] = maskIdNumber




    ### Extract the end state prior to a SN (isolated, only SMT, CEE, and merger) 

    # Get Merger data from SPs
    isMerger = SPs['Merger'][()] == 1
    ends_preSN = ((st1_f < 13) | (st1_f == 16)) & ((st2_f < 13) | (st2_f == 16))  # 16 for CHE which are MS stars for this purpose
    maskSp_MergesPreSn = isMerger & ends_preSN # don't include post SN mergers
    # Get CEE data from CEs
    maskCePreSN = ((st1_ce < 13) | (st1_ce == 16)) & ((st2_ce < 13) | (st2_ce == 16))
    sdsCePreSN = ceSeeds[maskCePreSN]
    maskSp_CePreSN = np.in1d(spSeeds, sdsCePreSN)
    # Combine mergers and CEE
    maskSp_MergerOrCeePreSn = maskSp_MergesPreSn | maskSp_CePreSN
    maskSp_CePreSN_nonMerging = maskSp_CePreSN & ~maskSp_MergesPreSn

    # Get interacting non-mergers by grabbing all interacters preSN, and removing any from the merger mask above
    maskMt_InteractionsPreSn = ((st1_mt < 13) | (st1_mt == 16)) & ((st2_mt < 13) | (st2_mt == 16)) # all MT events that occur prior to a SN
    # all MT events that occur prior to a SN
    sdsMt_InteractsPreSn = np.unique(mtSeeds[maskMt_InteractionsPreSn]) # the seeds of everything that interacts prior to a SN
    maskSp_InteractPreSn = np.in1d(spSeeds, sdsMt_InteractsPreSn)
    maskSp_InteractOnlyStable = maskSp_InteractPreSn & ~maskSp_MergerOrCeePreSn

    # Everything else is a non-interactor - combine the masks together
    end_state_mask_list = [maskSp_MergesPreSn, maskSp_CePreSN_nonMerging, maskSp_InteractOnlyStable] 
    checkMasksMutuallyExclusive(end_state_mask_list) # check that things look good
    end_state_masks = np.zeros_like(spSeeds)
    for ii, mask in enumerate(end_state_mask_list):
        maskIdNumber = ii + 1
        end_state_masks[mask] = maskIdNumber




    ### Define the output categories 

    # Shape based on CEE/Merger, interacting, and isolated

    state_labels = ['Isolated', 'Merger', 'CEE', 'Only SMT']
    state_markers = ['x', '*', 'o', '^']
    state_colors = ['black', 'red', 'blue', 'green']

    # Color based on stellar type of 1st interaction - black for no interaction - color_labels are just the stellar type integers
    stype_colors = ['black', 'darkorange', 'red', 'blue', 'green', 'chocolate', 'fuchsia']

    # Combine related array types for the output
    params = [q, rP, m1, mT]
    stype_masks_colors = [stype_masks, stype_colors]
    state_masks_markers_labels_colors = [mt1_state_masks, end_state_masks, state_labels, state_markers, state_colors]

    return params, stype_masks_colors, state_masks_markers_labels_colors


