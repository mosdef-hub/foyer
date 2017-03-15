import mbuild as mb

compound = mb.load('1-octanol.mol2')
compound.save('1-octanol-mb.mol2', forcefield_name='trappe-ua')
