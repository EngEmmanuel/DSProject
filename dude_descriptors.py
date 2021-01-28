# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 20:37:19 2021

@author: webma
"""

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
import pandas as pd
import gzip



dude_target_dict = {
    'CAH2': 'CHEMBL205',
    'JAK2': 'CHEMBL2971',
    'ACES': 'CHEMBL220',
    'VGFR2': 'CHEMBL4142',
    'BACE1': 'CHEMBL4822',
    'EGFR': 'CHEMBL203',
    'AA2AR': 'CHEMBL251',
    'FA10': 'CHEMBL244'
}

inf = gzip.open('data/cah2/actives_final.sdf.gz')
gzactives = Chem.ForwardSDMolSupplier(inf)
actives = [x for x in gzactives if x is not None]
#for mol in actives:
#   print(mol.GetNumAtoms())

inf = gzip.open('data/cah2/decoys_final.sdf.gz')
gzdecoys = Chem.ForwardSDMolSupplier(inf)
decoys = [x for x in gzdecoys if x is not None]

#for mol in decoys:
#   print(mol.GetNumAtoms())

descriptors_actives = {}

descriptors_actives['n_atoms'] = {'dude_cah2_a' + str(i): Descriptors.HeavyAtomCount(actives[i]) for i in range(0,len(actives))}
descriptors_actives['rotatable_bonds'] = {'dude_cah2_a' + str(i): Descriptors.NumRotatableBonds(actives[i]) for i in range(0,len(actives))}
descriptors_actives['logp'] = {'dude_cah2_a' + str(i): Descriptors.MolLogP(actives[i]) for i in range(0,len(actives))}
descriptors_actives['molecular_weight'] = {'dude_cah2_a' + str(i): Descriptors.ExactMolWt(actives[i]) for i in range(0,len(actives))}
descriptors_actives['hb_donors'] = {'dude_cah2_a' + str(i): Descriptors.NumHDonors(actives[i]) for i in range(0,len(actives))}
descriptors_actives['hb_acceptors'] = {'dude_cah2_a' + str(i): Descriptors.NumHAcceptors(actives[i]) for i in range(0,len(actives))}

descriptors_actives = pd.DataFrame.from_dict(descriptors_actives)

descriptors_decoys = {}

descriptors_decoys['n_atoms'] = {'dude_cah2_d' + str(i): Descriptors.HeavyAtomCount(decoys[i]) for i in range(0,len(decoys))}
descriptors_decoys['rotatable_bonds'] = {'dude_cah2_d' + str(i): Descriptors.NumRotatableBonds(decoys[i]) for i in range(0,len(decoys))}
descriptors_decoys['logp'] = {'dude_cah2_d' + str(i): Descriptors.MolLogP(decoys[i]) for i in range(0,len(decoys))}
descriptors_decoys['molecular_weight'] = {'dude_cah2_d' + str(i): Descriptors.ExactMolWt(decoys[i]) for i in range(0,len(decoys))}
descriptors_decoys['hb_donors'] = {'dude_cah2_d' + str(i): Descriptors.NumHDonors(decoys[i]) for i in range(0,len(decoys))}
descriptors_decoys['hb_acceptors'] = {'dude_cah2_d' + str(i): Descriptors.NumHAcceptors(decoys[i]) for i in range(0,len(decoys))}

descriptors_decoys = pd.DataFrame.from_dict(descriptors_decoys)
