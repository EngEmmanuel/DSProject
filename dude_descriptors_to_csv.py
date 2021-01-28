# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 20:37:19 2021

@author: webma
"""

from rdkit import Chem

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

import pandas as pd
import gzip
import sqlite3
import oddt
import os
from oddt import virtualscreening
from oddt.virtualscreening import electroshape



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

concatenated = pd.DataFrame()
for target in dude_target_dict:

    inf = gzip.open('data/' + target + '/actives_final.sdf.gz')
    gzactives = Chem.ForwardSDMolSupplier(inf)
    actives = [x for x in gzactives if x is not None]
    #for mol in actives:
    #   print(mol.GetNumAtoms())

    dfa = pd.DataFrame()
    dfa['target_chembl_id'] = [ dude_target_dict[target] for _ in actives ]
    dfa['dude_name'] = [ target for _ in actives ]
    dfa['active'] = 'active'
    dfa['uniquekey'] = ['dude_cah2_a' + str(i) for i in range(0,len(actives))]
    dfa['RDKit_Molecule'] = actives
    
    inf = gzip.open('data/' + target + '/decoys_final.sdf.gz')
    gzdecoys = Chem.ForwardSDMolSupplier(inf)
    decoys = [x for x in gzdecoys if x is not None]

    dfd = pd.DataFrame()
    dfd['target_chembl_id'] = [ dude_target_dict[target] for _ in decoys ]
    dfd['dude_name'] = [ target for _ in decoys ]
    dfd['active'] = 'decoy'
    dfd['uniquekey'] = ['dude_cah2_d' + str(i) for i in range(0,len(decoys))]
    dfd['RDKit_Molecule'] = decoys
    
    concatenated = pd.concat([dfa, dfd])
    
    #for mol in decoys:
    #   print(mol.GetNumAtoms())
    


# Helper function to compute descriptors for a single molecule – from Fergus B
def compute_descriptors(molecule):
    descriptors = {d[0]: d[1](molecule) for d in Descriptors.descList}
    descriptors = pd.Series(descriptors)
    return descriptors

descriptors = concatenated['RDKit_Molecule'].apply(compute_descriptors)

df_w_rdkit_desc = pd.concat([concatenated, descriptors], axis=1)
df_w_rdkit_desc.to_csv('data/10_dude_data_plus_rdkit_descriptors.csv', index=False)


# descriptors_actives = {}

# descriptors_actives['n_atoms'] = {'dude_cah2_a' + str(i): Descriptors.HeavyAtomCount(actives[i]) for i in range(0,len(actives))}
# descriptors_actives['rotatable_bonds'] = {'dude_cah2_a' + str(i): Descriptors.NumRotatableBonds(actives[i]) for i in range(0,len(actives))}
# descriptors_actives['logp'] = {'dude_cah2_a' + str(i): Descriptors.MolLogP(actives[i]) for i in range(0,len(actives))}
# descriptors_actives['molecular_weight'] = {'dude_cah2_a' + str(i): Descriptors.ExactMolWt(actives[i]) for i in range(0,len(actives))}
# descriptors_actives['hb_donors'] = {'dude_cah2_a' + str(i): Descriptors.NumHDonors(actives[i]) for i in range(0,len(actives))}
# descriptors_actives['hb_acceptors'] = {'dude_cah2_a' + str(i): Descriptors.NumHAcceptors(actives[i]) for i in range(0,len(actives))}

# descriptors_actives = pd.DataFrame.from_dict(descriptors_actives)

# descriptors_decoys = {}

# descriptors_decoys['n_atoms'] = {'dude_cah2_d' + str(i): Descriptors.HeavyAtomCount(decoys[i]) for i in range(0,len(decoys))}
# descriptors_decoys['rotatable_bonds'] = {'dude_cah2_d' + str(i): Descriptors.NumRotatableBonds(decoys[i]) for i in range(0,len(decoys))}
# descriptors_decoys['logp'] = {'dude_cah2_d' + str(i): Descriptors.MolLogP(decoys[i]) for i in range(0,len(decoys))}
# descriptors_decoys['molecular_weight'] = {'dude_cah2_d' + str(i): Descriptors.ExactMolWt(decoys[i]) for i in range(0,len(decoys))}
# descriptors_decoys['hb_donors'] = {'dude_cah2_d' + str(i): Descriptors.NumHDonors(decoys[i]) for i in range(0,len(decoys))}
# descriptors_decoys['hb_acceptors'] = {'dude_cah2_d' + str(i): Descriptors.NumHAcceptors(decoys[i]) for i in range(0,len(decoys))}

# descriptors_decoys = pd.DataFrame.from_dict(descriptors_decoys)

# chembldb.execute(query)

# query = """INSERT INTO dude_chembl (dude_name, chembl_id) VALUES
# ('CAH2', 'CHEMBL205'),
# ('JAK2', 'CHEMBL2971'),
# ('ACES', 'CHEMBL220'),
# ('VGFR2', 'CHEMBL4142'),
# ('BACE1', 'CHEMBL4822'),
# ('EGFR', 'CHEMBL203'),
# ('AA2AR', 'CHEMBL251'),
# ('FA10', 'CHEMBL244')
# """
# chembldb.execute(query)
# chembldb.execute("COMMIT")


