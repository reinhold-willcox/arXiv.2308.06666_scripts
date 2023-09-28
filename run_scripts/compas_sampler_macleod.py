import os, sys
import numpy as np
import scipy.stats as ss
import subprocess

### Set run options 
run_on_hpc = True # Run on slurm based cluster HPC # keep set to True on the tracked version
if run_on_hpc:
    nCores = 25
    print("Running on HPC with {} cores".format(nCores))
else:
    nCores = 12 # Set to the number of local cores to use
    print("Running locally with {} cores".format(nCores))
logfile_defs = 'logfile_defs.txt' #None

### Set sampling parameters
nSystems = 1e5 # Decide later whether / how to split into singles or binaries
seed_base = 0
useBinFrac = False 
binFrac = 0.5 # only applies if useBinFrac is True


primaryOutputDir = 'masterOutput/'

###############################################################################
# ##
# ## Define the model and model parameters
# ##
# ######################################################################################

### For model variations that come from compas input parameters, set those as dictionaries below
modelConfigMap = {

        #'fgamma0.00_betaComp':
        #     { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #       '--mass-transfer-jloss-macleod-linear-fraction': 0.0 },
        #'fgamma0.25_betaComp':
        #     { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #       '--mass-transfer-jloss-macleod-linear-fraction': 0.25 },
        #'fgamma0.50_betaComp':
        #     { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #       '--mass-transfer-jloss-macleod-linear-fraction': 0.5 },
        #'fgamma0.75_betaComp':
        #     { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #       '--mass-transfer-jloss-macleod-linear-fraction': 0.75 },
        #'fgamma1.00_betaComp':
        #    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #      '--mass-transfer-jloss-macleod-linear-fraction': 1.0 },
        #'fgamma1.00_beta0.00':
        #    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #      '--mass-transfer-jloss-macleod-linear-fraction': 1.0 ,
        #      '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 0.0 },
        #'fgamma0.50_beta0.00':
        #    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
        #      '--mass-transfer-jloss-macleod-linear-fraction': 0.5 ,
        #      '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 0.0 },
        #'beta0.00':
        #    { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 0.0 },
        #'beta0.25':
        #    { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 0.25 },
        #'beta0.50':
        #    { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 0.5 },
        #'beta0.75':
        #    { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 0.75 },
        #'beta1.00':
        #    { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
        #      '--mass-transfer-fa' : 1.0 },
        #'zetaMsInf':
        #    { '--zeta-main-sequence': 10000 },
        #'zetaGiantInf':
        #    { '--zeta-radiative-envelope-giant': 10000 },
        #'Tconv3.73': # T = 5730
        #    {'--envelope-state-prescription' : 'FIXED_TEMPERATURE'},
        #'Tconv3.65': 
        #    {'--envelope-state-prescription' : 'FIXED_TEMPERATURE'
        #        '--convective-envelope-temperature-threshold': 4470},
        #'Tconv3.90':
        #    {'--envelope-state-prescription' : 'FIXED_TEMPERATURE'
        #        '--convective-envelope-temperature-threshold': 7940},
        #'Tconv3.60':
        #    {'--envelope-state-prescription' : 'FIXED_TEMPERATURE'
        #        '--convective-envelope-temperature-threshold': 3980},
        'qCritClaeys':
            {'--critical-mass-ratio-prescription' : 'CLAEYS'},
        'qCritGe20':
            {'--critical-mass-ratio-prescription' : 'GE20'},
        'qCritGe20IC':
            {'--critical-mass-ratio-prescription' : 'GE20_IC'},

}

# Sample in different T convs
#for jj, Tconv in enumerate( ['Tconv3.73', 'Tconv3.65']):
#    moddict = dict()
#    modname = 'zetaDefault_{}'.format(Tconv)
#    moddict.update({'--envelope-state-prescription' : 'FIXED_TEMPERATURE'})
#    if Tconv == 'Tconv3.65':
#        moddict.update({'--convective-envelope-temperature-threshold': 4470})
#    modelConfigMap.update({modname: moddict.copy()})



# Sample extra thermal C's
#primaryOutputDir = 'extra_thermal_Cs/'
#for jj, zeta in enumerate(['Default', 'Ge20']): #, 'Ge20IC']):
#    for cc, beta in enumerate(['CompC03', 'CompC30']): # ['0.0', '0.5', '1.0', 'Comp']:
#
#        moddict = dict()
#        modname = 'zeta{}_beta{}'.format(zeta, beta)
#
#        if jj > 0:
#            zeta_arg = ['', 'GE20', 'GE20_IC'][jj]
#            moddict.update({'--stellar-zeta-prescription' : zeta_arg})
#
#        if beta not in ['Comp', 'CompC03', 'CompC30']:
#            moddict.update(
#              { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
#                '--mass-transfer-fa' : float(beta) },
#              )
#        elif beta in ['CompC03', 'CompC30']:
#            Cval = [3, 30][cc]
#            moddict.update( { '--mass-transfer-thermal-limit-C' : Cval })
#
#        if beta != '1.0':
#            for fg in ['0.00', '0.25', '0.50', '0.75', '1.00']: # probably don't even need all these 
#                moddict.update( 
#                    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
#                      '--mass-transfer-jloss-macleod-linear-fraction': float(fg) }
#                    )
#                newmodname = modname+'_fgamma{}'.format(fg)
#                modelConfigMap.update({newmodname: moddict.copy()})
#        else:
#            modelConfigMap.update({modname: moddict.copy()})




# Sample coarsely in f_gamma and beta
for jj, zeta in enumerate(['Default']):
    for beta in ['0.0', '0.5', '1.0', 'Comp']:

        moddict = dict()
        modname = 'zeta{}_beta{}'.format(zeta, beta)

        if jj > 0:
            zeta_arg = ['', 'GE20', 'GE20_IC'][jj]
            moddict.update({'--stellar-zeta-prescription' : zeta_arg})

        if beta != 'Comp':
            moddict.update(
              { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
                '--mass-transfer-fa' : float(beta) },
              )
        if beta != '1.0':
            for fg in ['0.00', '0.25', '0.50', '0.75', '1.00']: # probably don't even need all these 
                moddict.update( 
                    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
                      '--mass-transfer-jloss-macleod-linear-fraction': float(fg) }
                    )
                newmodname = modname+'_fgamma{}'.format(fg)
                modelConfigMap.update({newmodname: moddict.copy()})
        else:
            modelConfigMap.update({modname: moddict.copy()})


## Sample finely in f_gamma and beta
#primaryOutputDir = 'friviolous_extra_gammas/'
#for jj, zeta in enumerate(['Default', 'Ge20']):
#    for beta in ['0.0', '0.5', '1.0', 'Comp']:
#
#        moddict = dict()
#        modname = 'zeta{}_beta{}'.format(zeta, beta)
#
#        if jj > 0:
#            zeta_arg = ['', 'GE20', 'GE20_IC'][jj]
#            moddict.update({'--stellar-zeta-prescription' : zeta_arg})
#
#        if beta != 'Comp':
#            moddict.update(
#              { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
#                '--mass-transfer-fa' : float(beta) },
#              )
#        if beta != '1.0':
#            for fg in ['0.00', '0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90', '1.00']:
#                moddict.update( 
#                    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
#                      '--mass-transfer-jloss-macleod-linear-fraction': float(fg) }
#                    )
#                newmodname = modname+'_fgamma{}'.format(fg)
#                modelConfigMap.update({newmodname: moddict.copy()})
#        else:
#            modelConfigMap.update({modname: moddict.copy()})



## Sample in different zetas, across some of the betas
#for jj, zeta in enumerate(['Ge20', 'Ge20IC']):
#    for beta in ['0.0', '0.5', '1.0', 'Comp']:
#        moddict = dict()
#        modname = 'zeta{}_beta{}'.format(zeta, beta)
#
#        zeta_arg = ['GE20', 'GE20_IC'][jj]
#        moddict.update({'--stellar-zeta-prescription' : zeta_arg})
#
#        if beta != 'Comp':
#            moddict.update(
#              { '--mass-transfer-accretion-efficiency-prescription' : 'FIXED',
#                '--mass-transfer-fa' : float(beta) },
#              )
#        if beta != '1.0':
#            for fg in [ '0.00', '0.25', '0.50', '0.75', '1.00']:
#                
#                moddict.update(
#                    { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
#                      '--mass-transfer-jloss-macleod-linear-fraction': float(fg) }
#                    )
#                newmodname = modname+'_fgamma{}'.format(fg)
#                modelConfigMap.update({newmodname: moddict.copy()})
#        else:
#            modelConfigMap.update({modname: moddict.copy()})




#for key, val in modelConfigMap.items():
#    print(key)
#    for k, v in val.items():
#        print(k, v)
#    print()
#print(len(modelConfigMap))
#print(done)

# generate pair-wise gamma, max time scrips
#for tim in [1, 3, 10, 30, 100, 300]:
#    for gam in ['0.0', '0.5', '1.0']:
#
#        modname = 'tmax{:03}_gamma{}'.format(tim, gam)
#        moddict = \
#            { '--mass-transfer-angular-momentum-loss-prescription' : 'MACLEOD_LINEAR',
#              '--mass-transfer-jloss-macleod-linear-fraction': float(gam),
#              '--maximum-evolution-time': float(tim),
#              }
#        modelConfigMap.update({modname: moddict})


index = int(sys.argv[1]) # numerical arguments following the python script call decide which model to use
nMods = len(list(modelConfigMap.keys()))
model = list(modelConfigMap)[index % nMods]

if index < nMods:
    should_create_new_grid_file = False # use moe di stefano grid
    moddir = 'moedistef' 
else:
    should_create_new_grid_file = True # run with moe and distef for first half
    moddir = 'uniformIC' 

existing_gridfile = 'grid_moedistefano_nSamples100000.txt' # Fill in grid name

print('Using model: {}'.format(model))


### Set paths to the exe dir and output folder
#compas_root = '/home/rwillcox/astro/compas/COMPAS_MACLEOD/' # Ensure this is set to the correct dir
if run_on_hpc:
    compas_root = '/home/rwillcox/astro/compas/COMPAS_MACLEOD/' # Ensure this is set to the correct dir
    output_path = '/fred/oz101/rwillcox/macleodgamma/{}/{}/{}/'.format(primaryOutputDir, moddir, model) 
else:
    compas_root = '/home/rwillcox/astro/compas/COMPAS/' # Ensure this is set to the correct dir
    output_path= '/home/rwillcox/astro/MacLeod_gamma/data/{}'.format(model)


# Set non-sampled parameters, which do not vary between runs, here
command_line_args = { 
    '--evolve-unbound-systems' : 'TRUE',
    #'--use-critical-mass-ratio-HG' : 'TRUE',
    #'--use-critical-mass-ratio-MS-high-mass' : 'TRUE',
    #'--use-critical-mass-ratio-giant' : 'TRUE',
}

# Update the command line args as needed
# logfile defs not needed now that mergers have been added completely
#if logfile_defs is not None: 
#    command_line_args.update({'--logfile-definitions': logfile_defs}) 
command_line_args.update(modelConfigMap[model])   # Update with the specific parameters for this run



# Set number of singles and binaries
if useBinFrac: # Can set the number of single stars directly from the binary fraction
    nBinaries = int(np.floor(nSystems*binFrac))  
    nSingles = int(np.floor(nSystems*(1-binFrac)))  
else:
    nBinaries = nSystems
    nSingles = 0


# Set sampled parameters here
def sample_parameters(model=None, isBinary=True):
    """
    This function should return a single sample for each parameter,
    and is expected to be called potentially many times.
    Check later if slow.
    """

    # KroupaIMF, with alpha=-2.3
    m1Min = 5
    if model == 'm1Min5':
        m1Min = 5
    m1Max = 100
    slope = -2.3
    m1 = sample_power_law(slope=slope, xmin=m1Min, xmax=m1Max)

    if isBinary:
        ### Secondary mass sampled from mass ratio
        # default values
        qMin = 0.1
        qMax = 1.0
        qSlope = 0.0
        if model == 'qMin0.5':
            qMin = 0.5
        if model == 'qMax0.5':
            qMax = 0.5
        if model == 'qLow':
            qSlope = -1
        if model == 'qHigh':
            qSlope = 1
        q = sample_power_law(slope=qSlope, xmin=qMin, xmax=qMax)
        m2 = m1*q

        if (model == 'zapartas19') or (model == 'zapartasPZams'):
            # Here, sample in period
            logPMin = 0.15
            logPMax = 3.5
            slope = -0.55 if (m1>15) else 0  #slope = 0 #-0.55 if (m1>15) else 0
            logP = sample_power_law(slope=slope, xmin=logPMin, xmax=logPMax)
            P = np.power(10, logP) # days
            Pyr = P/365.25
            a = np.cbrt( (m1+m2) * Pyr*Pyr) # AU
        else: 
            ### Sample semi-major axis
            # default values
            logAMin = -2
            logAMax = 2
            logASlope = 0
            if model == 'aMin10AU':
                logAMin = 1
            if model == 'aMax10AU':
                logAMax = 1
            if model == 'aLow':
                logASlope = -1
            if model == 'aHigh':
                logASlope = 1
            logA = sample_power_law(slope=logASlope, xmin=logAMin, xmax=logAMax)
            a = np.power(10, logA)

        sampled_params = {
            '--initial-mass-1' : m1,
            '--initial-mass-2' : m2,
            '-a' : a,
        }

    else: # Just "single" stars - make the secondary small and separation large
        m2 = 1.0
        a = 1000 # AU

        sampled_params = {
            '--initial-mass-1' : m1,
            '--initial-mass-2' : m2,
            '-a' : a
        }
    
    return sampled_params




###############################################################################
# ##
# ## Leave the lower parts untouched, only set configurations above 
# ##
# ######################################################################################

def main():
    # Create directory structure
    create_dir_structure(output_path)

    # Print the model name and set the random seed
    np.random.seed(seed_base) 

    # Iterate over the number of cores
    nBinariesPerCore = int(np.ceil(nBinaries/nCores)) # ceiling so you have >= the number of requested binaries
    nSinglesPerCore = int(np.ceil(nSingles/nCores)) # ceiling so you have >= the number of requested binaries

    # Create the file structure, run compas, and recombine the outputs when the runs complete
    pIDs = createFileStructureAndStartRuns(nBinariesPerCore, nSinglesPerCore)
    recombineOutputs(pIDs)

def createFileStructureAndStartRuns(nBinariesPerCore, nSinglesPerCore):
    # Create all the grid files
    process_ids = []

    for index in range(nCores):

        container = 'batch_{}'.format(index)
        gridfile = os.path.join(output_path, 'grid_{}.csv'.format(index, index))
        logfile = os.path.join(output_path, 'logs/log_{}.txt'.format(index))
        compas_exe = os.path.join(compas_root, 'src/COMPAS') # Location of COMPAS executable 

        # Sample the correct parameters and put into grid file
        if should_create_new_grid_file:
            create_new_grid_file(gridfile, nBinariesPerCore, nSinglesPerCore)
        else:
            extract_from_existing_grid_file(existing_gridfile, gridfile, nCores, index)

        # Run COMPAS on the appropriate grid
        slurm_file = os.path.join(output_path, 'slurms/slurm_{}.sh'.format(index))
        nSystemsPerCore = nBinariesPerCore + nSinglesPerCore
        create_slurm_file(slurm_file, compas_exe, output_path, gridfile, container, logfile, nSystemsPerCore, index, options=command_line_args)
        
        if run_on_hpc:
            process = subprocess.Popen("sbatch " + slurm_file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # run command in the background
            stdout = process.communicate()[0].decode("utf-8")
            pid = stdout.split()[3]
            process_ids.append(pid)
        else:
            subprocess.Popen("chmod +x " + slurm_file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # run command in the background
            process = subprocess.Popen("/bin/sh " + slurm_file + ' &', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # run command in the background
            process_ids.append(process)
        
    return process_ids




def recombineOutputs(process_ids):

    # Check that COMPAS runs have completed
    subprocessesCompleted = False
    nPIDs = len(process_ids)
    frac_completed = 0
    if run_on_hpc:
        while not subprocessesCompleted:
            # Get list of ongoing process IDs
            proc = subprocess.Popen("squeue -u rwillcox | awk '{print $1}'", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # run command in the background
            ongoing_processes = proc.communicate()[0].decode("utf-8").split()
            nIDsRemaining = sum(p in ongoing_processes for p in process_ids) # If none of the ongoing process was started above:
            frac_done = (nPIDs - nIDsRemaining) / nPIDs
            if frac_done > frac_completed:
                print("Completed {:.1f}%".format(frac_done*100))
                frac_completed = frac_done 
            if nIDsRemaining == 0:
                subprocessesCompleted = True
    else:
        # Process_ids are just a list of the background processes, need to run communicate() now
        for pp, process in enumerate(process_ids):
            # Running communicate waits for the process to finish, so this might be enough
            process.communicate()[0]
            frac_done = (pp+1)/nPIDs
            print("Completed {:.1f}%".format(frac_done*100))

    # Recombine them into one output hdf5 file
    h5copyFile = os.path.join(compas_root, 'utils/h5copy.py') # Location of COMPAS h5copy File
    shell_cmd = "python3 {h5copyFile} -o {output_path}/COMPAS_Output.h5 {output_path}/*".format(h5copyFile=h5copyFile, output_path=output_path)
    subprocess.run(shell_cmd, shell=True)



# Power law sampler 
def sample_power_law(slope, xmin, xmax, size=1):
    """
    Need to define my custom power law sampler 
    since Scipy can't handle negative exponents...

    Add functionality for size later if slow
    """
    ### Do inverse CDF sampling of power law
    # If slope is -1, there are division errors
    if slope == -1:
        slope = -.99
    # If slope is neg and some of [xmin, xmax] < 0, this will create a nan. Keep the slope but translate the x values up to positive values, then subtract off afterward
    if slope < 0:
        offset = abs(xmin) - xmin # will be 0 for non-negative xmin
        xmin += offset
        xmax += offset
    else:
        offset = 0
    
    return np.power(np.random.uniform(0, 1, size) * (np.power(xmax, 1+slope) - np.power(xmin, 1+slope)) + np.power(xmin, 1+slope), 1/(1+slope))[0] -offset




def create_dir_structure(outpath):
    if os.path.exists(outpath):
        raise FileExistsError("Output path exists, aborting now\n{}".format(outpath))
    else:
        os.makedirs(outpath)
        # make slurm and log subdirectories
        os.makedirs(os.path.join(outpath, 'slurms'))
        os.makedirs(os.path.join(outpath, 'logs'))

def create_slurm_file(slurm_file, exe, outpath, gridfile, container, logfile, nSystemsPerCore, index, options=dict()):
    """ I'm calling any file containing the run line a slurm file, even if running on non-hpc where slurm is irrelevant"""
    if run_on_hpc:
        total_seconds = nSystemsPerCore * .3 # roughly 0.1s/system
        hh = int(total_seconds // 3600)
        mm = int((total_seconds - hh*3600 )//60)
        ss = int(round(total_seconds - hh*3600 - mm*60))
        slurm_string_base = "#!/bin/bash\n#SBATCH --mem-per-cpu=1024\n#SBATCH --output=slurm_output.out\n#SBATCH --time=0-{hh}:{mm}:{ss}\n".format(hh=hh,mm=mm,ss=ss)
    else:
        slurm_string_base = ""
    seed = seed_base + index*nSystemsPerCore
    options_string = " ".join(["{} {}".format(key, val) for key, val in options.items()])
    slurm_string = "{base}{exe} --grid {grid} --output-container {container} --output-path {outpath} --random-seed {seed} {options} > {logfile}".format(
            base=slurm_string_base, exe=exe, grid=gridfile, container=container, outpath=outpath, seed=seed, options=options_string, logfile=logfile)
    with open(slurm_file, 'w') as fwrite:
        fwrite.write(slurm_string)

def create_new_grid_file(gridfile, nBinaries, nSingles):
    with open(gridfile, 'w') as fwrite:
        for jj in range(nBinaries):
            fwrite.write(" ".join(["{} {}".format(key, val) for key, val in sample_parameters(model=model).items()])+"\n")
        for kk in range(nSingles):
            fwrite.write(" ".join(["{} {}".format(key, val) for key, val in sample_parameters(model=model, isBinary=False).items()])+"\n")


def extract_from_existing_grid_file(existing_gridfile, new_gridfile, nCores, index):
    
    # Extract every nth line of existing_gridfile into new_gridfile, where n is the index
    # this avoids having to iterate over the input lines twice just to get the total line number

    with open(existing_gridfile) as fread:
        with open(new_gridfile, 'w') as fwrite:
            for line_num, line in enumerate(fread):
                if (line_num % nCores) == index:
                    fwrite.write(line)


def test():
    with open('testing/testout_{}_{}.txt'.format(outdir, index), 'w') as fwrite:
        [fwrite.write("{} : {}\n".format(key,val)) for key, val in command_line_args.items()]
        fwrite.write("\n")
        for jj in range(10):
            fwrite.write(" ".join(["{} {}".format(key, val) for key, val in sample_parameters(model=model).items()])+"\n")

if __name__ == "__main__":
    main()
