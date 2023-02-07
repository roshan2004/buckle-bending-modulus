### Author: Jackson Crowley
#
### Date: 07/02/2023
#
### Version: 1.0
#
# Description: After following the protocol for the calculation of bending modulus via the buckling
# method, outlined on our team GitHub page (https://github.com/MMSB-MOBI/buckle-bending-modulus), the user will have performed two production
# simulatons: one for the calculation of length, where the bilayer is free to move only in x,
# and a second where the biayer is buckled in x, exerting force on the "wall" of the box in x.
#
# Usage: This program requires two flags:
#               -l being the gmx energy file (.edr) output from the length simulation
#               -f being the gmx energy file (.edr) output from the force/buckled simulation
#
#        The optional flag -v will write out a more human readable output, whereas not signalling
#        will output a simple list of the full Kc and Kc Error.
#
#        The optional flag -k will keep the outputs from the "gmx energy" call.
                                   

import os
import subprocess
import argparse
import numpy as np
import math

parser = argparse.ArgumentParser(description = "Calculate the Bending Modulus, Kc, from a buckled membrane. \
                                 requires two inputs, -l being the GROMACS energy output from the system used \
                                to calculate length in the x-dimension, and -f being the energy output corresponding \
                                 to the system buckled in the same dimension.",

epilog = "Note that this script is simply an implementation of the equation derived in Hu, 2013 (doi: 10.1063/1.4808077).")

# -c INPUTFILE -o OUTPUTFILE -x DIM -y DIM -z DIM --help
parser.add_argument("-l", "--length", help = "<.edr> GROMACS energy file from system length simulation.", required = True)
parser.add_argument("-f", "--force", help = "<.edr> GROMACS energy file from force exerted by buckle", required = True)
parser.add_argument("-v", "--verbose", default = False, action='store_true', help = "if indicated, -v will output the result in a more human readable form.")
parser.add_argument("-k", "--keep", default = False, action='store_true', help = "Will keep the created GROMACS energy output files, if flagged, which are deleted otherwise.")
args = parser.parse_args()

def calculate_bending(L, L_err, Lx, Ly, Lz, P_xx, P_xx_error, T):
    # unit corrections
    L = L*1e-9
    L_err = L_err*1e-9
    Lx = Lx*1e-9
    Ly = Ly*1e-9
    Lz = Lz*1e-9
    P_xx = P_xx*1e5
    P_xx_error = P_xx_error*1e5

    # constants
    kb = 1.3806488e-23
    bi_constant  = [1 ,0.5,0.28125,0.1640625,0.0970458984375,0.05767822265625,0.034286499023438,0.020331859588623,0.012007629498839,0.007054503075778,0.004119324556086]
    di_constant = [ 1,0.625,0.421875,0.2880859375,0.1959228515625,0.13201904296875,0.087979316711426,0.057950466871262,0.037723844870925,0.024270858848468,0.01543510783813]

    # calculate the b coefficient from strain and correction
    strain = (L-Lx)/L
    
    sum_b_strain = 0
    for i in range( 0 , len(bi_constant)):
        sum_b_strain +=  bi_constant[i]*math.pow(strain, i) 
        
    coef_bi = L*L / ( 4 * math.pi*math.pi *Ly*sum_b_strain )

    # force and correction
    Fx= (Ly*Lz)*P_xx 

    sum_d_strain = 0
    for i in range( 0 , len(di_constant)):
        sum_d_strain +=  di_constant[i]*math.pow(strain, i) 

    coef_term = -(3*kb*T)/(2*L)    

    Fx_corrected = Fx + ( coef_term * sum_d_strain)

    # calculate bending modulus from force and b coefficient
    kappa_j = coef_bi * Fx_corrected
    bending_modulus =  kappa_j/(kb*T)

    # error propogation and calculation
    coef_error = coef_bi*Ly*Lz/(L*L)
    errL2 = math.sqrt( 2)*L*L_err
    errL2Pxx = math.sqrt( (P_xx *P_xx * errL2 * errL2) + (L*L * L*L * P_xx_error * P_xx_error) )
    error_propagation = coef_error * errL2Pxx
    bending_modulus_err = error_propagation / ( T * kb)
    
    return(bending_modulus, bending_modulus_err)

## Read input files and extract data for calculation

subprocess.call('echo -e "Box-X\n0" | gmx energy -f ' + args.length + ' -o lx.xvg', shell=True, stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)
Lx_xvg = open('lx.xvg', 'r')
L_list = []
for line in Lx_xvg.readlines():
    if line.startswith(("#", "@")):
        continue
    else:
        linesplit = line.split()
        L_list.append(float(linesplit[1]))
        
L = np.average(L_list)
L_err = np.std(L_list)

subprocess.call('echo -e "Temperature\nBox-X\nBox-Y\nBox-Z\nPres-XX\n0"| gmx energy -f ' + args.force + ' -o force.xvg', shell=True, stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)

bending_xvg = open('force.xvg', 'r')
T_list = []
X_list = []
Y_list = []
Z_list = []
P_xx_list = []
for line in bending_xvg.readlines():
    if line.startswith(("#", "@")):
        continue
    else:
        linesplit = line.split()
        T_list.append(float(linesplit[1]))
        X_list.append(float(linesplit[2]))
        Y_list.append(float(linesplit[3]))
        Z_list.append(float(linesplit[4]))
        P_xx_list.append(float(linesplit[5]))
        
T = np.average(T_list)
Lx = np.average(X_list)
Ly = np.average(Y_list)
Lz = np.average(Z_list)
P_xx = np.average(P_xx_list) - 1 # we subtract 1 due to the water pressure
P_xx_halfway = len(P_xx_list)/2
np.average(P_xx_list[:len(P_xx_list)//2])
P_xx_err = (np.abs(np.average(P_xx_list[:len(P_xx_list)//2]) - np.average(P_xx_list[(len(P_xx_list)//2):])))


bending_calc = calculate_bending(L, L_err, Lx, Ly, Lz, P_xx, P_xx_err, T)
bending_modulus = bending_calc[0]
bending_modulus_err = bending_calc[1]


if args.verbose == True:
    print('\nKc= {:.3f} ± {:.3f} kBT'.format(bending_modulus, bending_modulus_err))
else:
    print(bending_calc)


if args.keep == False:
    os.remove("lx.xvg")
    os.remove("force.xvg")
else:
    print("\n\tgmx energy outputs for length and force simulations saved as lx.xvg and force.xvg, respectively.\n")
