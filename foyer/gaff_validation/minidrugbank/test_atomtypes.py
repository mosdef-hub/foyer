# coding: utf-8
import foyer
import mbuild
import parmed as pmd
import glob
import re

# Define the forcefield
ff = foyer.Forcefield(forcefield_files='../../forcefields/gaff.xml')

# Find the test files

load_path = 'ante_types/'
save_path = 'foyer_types/'
log_path = 'compare_antechamber-foyer.log'
files = glob.glob(load_path+'/*.mol2')

with open(log_path,'w') as log:

    for filen in files:
        molname = re.sub('.mol2','',re.sub('^.*/','',filen))
        log.write('{}\n'.format(molname))
        print('Starting {} ... \n\n'.format(molname))
        # Load a molecule and extract antechamber types
        molecule = pmd.load_file(filen)
        types_ante = [ atom.type for atom in molecule ]
        # Convert to parmed structure, run foyer atomtyping 
        untyped = molecule.to_structure()
        # For now -- just check atom typing -- not the existence of parameters
        foyer_typed = ff.apply(untyped,assert_bond_params=False,assert_angle_params=False,assert_dihedral_params=False)
        # Regex to remove extra _ tags for 'extra' atomtypes
        types_foyer = [ re.sub('_.*$','',atom.type) for atom in foyer_typed ]
   
        # Compare
        assert len(types_foyer) == len(types_ante), 'Different numbers of atoms!'
        idxs = range(len(types_foyer))
    
        for idx,type_a,type_f in zip(idxs,types_ante,types_foyer):
            if type_a != type_f:
                log.write('idx: {:5d} Ante: {:2s}  Foyer: {:2s}\n'.format(idx,type_a,type_f))
 
