### Please note this script is not designed to work out of the box, but rather it is an overview
### of the workflow typically followed to create the desired systems for Kc calculation.
###
### DOPC will be used as a model lipid here, as well as in the provided .mdp files in /mdp

### 1. Build Bilayer with Insane.py
# x, y, and z size used to exceed possibility of finite size effects (doi: 10.1021/acs.jpcb.0c04253)
###############################################################################################
mkdir /1.init
python2 insane.py -l DOPC -sol W -o insane_DOPC.gro -p insane_DOPC.top -x 32.5 -y 8.5 -z 20

## Perform neutralisation (even if the lipids are uncharged, easier when looping over everything)
# insane_DOPC.top will need to be replaces with a topology appropriately referencing your .itps,
# I'll define a placeholder "topol.top"
gmx grompp -f mdp/em.mdp -c 1.init/insane_DOPC.gro -p topol.top -o neutralise.tpr -maxwarn 2
echo 3 | gmx genion -s neutralise.tpr -p ../topol.top -neutral -o neutral.gro

## Make index; note you may need to adjust mdps in presence of ions
echo -e "del 0-19\nr " + lipid + "\na W WF\nname 1 Solvent\na NA CL\nname 2 ION\nq " \
 | gmx make_ndx gmx make_ndx -f 1.init/neutral.gro

### 2. Energy Minimise
###############################################################################################
mkdir 2.em
gmx grompp -f mdp/em.mdp -c 1.init/neutral.gro -p topol.top -o 2.em/em.tpr -maxwarn 1
gmx mdrun -v -deffnm 2.em/em

### 3. Equilibration
###############################################################################################
mkdir 3.eq 
gmx grompp -f mdp/eq.mdp -c 2.em/em.gro -p topol.top -o 3.eq/eq.tpr
gmx mdrun -v -deffnm 3.eq/eq 

### 4. LENGTH PRODUCTION RUN (Lx)
# This is one of the two simulations required for eventual calculation.
###############################################################################################
mkdir 4.lx
gmx grompp -f mdp/lx.mdp -c 3.eq/eq.gro -p topol.top -o 4.lx/lx.tpr
gmx mdrun -v -deffnm 4.lx/lx 

### 5. Deform
# (note that this is dependent on step 3, not step 4)
###############################################################################################
mkdir 5.deform
gmx grompp -f mdp/deform.mdp -c 3.eq/eq.gro -p topol.top -o 5.deform/deform.tpr
gmx mdrun -v -deffnm 5.deform/deform

### 6. Buckle Equilibration
# now that the system has been buckled, we need to equilibrate it in it's restrained form
# before the production run.
###############################################################################################
mkdir 6.buckle_eq
gmx grompp -f mdp/buckle_eq.mdp -c 5.deform/deform.gro -p topol.top -o 6.buckle_eq/buckle_eq.tpr
gmx mdrun -v -deffnm 6.buckle_eq/buckle_eq

### 7. RESTRAINED FORCE PRODUCTION RUN (Fx)
# This is the other simulation required for Kc calculation.
###############################################################################################
mkdir 7.restrain
gmx grompp -f mdp/restrain.mdp -c 6.buckle_eq/buckle_eq.gro -p topol.top -o 7.restrain/restrain.tpr
gmx mdrun -v -deffnm 7.restrain/restrain

### Calculation:
# if steps 4 and 7 have completed, you can call (some variation of) this command to calculate Kc:
python calculate_kc.py -l 4.lx/lx.edr -f 7.restrain/restrain.edr 