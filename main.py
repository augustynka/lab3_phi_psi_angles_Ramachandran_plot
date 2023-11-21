from Bio.PDB import PDBList, PDBParser, DSSP, PPBuilder
from Bio.PDB.Polypeptide import is_aa
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os


def calculate_phi_psi_angles(structure, dssp_path):
    phi_psi_data = []
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp=dssp_path)
    for chain in model:
        polypeptides = PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            # pobieranie kątów phi i psi
            phi_psi_list = poly.get_phi_psi_list()
            for i, (phi, psi) in enumerate(phi_psi_list):
                residue = poly[i]
                if is_aa(residue.get_resname(), standard=True):
                    if phi and psi:
                        # pobieranie struktury drugorzędowej z DSSP
                        sec_structure = dssp[(chain.id, residue.id)][2]
                        phi_psi_data.append((np.degrees(phi), np.degrees(psi), sec_structure))
    return phi_psi_data


def plot_ramachandran(phi_psi_data):
    data = pd.DataFrame(phi_psi_data, columns=['phi', 'psi', 'structure'])
    palette = {'H': 'red', 'E': 'blue', 'C': 'green', 'T': 'orange', 'S': 'purple', 'G': 'pink', '-': 'gray'}
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='phi', y='psi', hue='structure', palette=palette, data=data)

    plt.xlim(-180, 180)
    plt.ylim(-180, 180)

    plt.xlabel('Phi (°)')
    plt.ylabel('Psi (°)')
    plt.title('Wykres Ramachandrana')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks([-180, -120, -60, 0, 60, 120, 180])
    plt.yticks([-180, -120, -60, 0, 60, 120, 180])
    plt.tight_layout()
    #plt.show()
    plt.savefig('wykres_ramachandrana.png')


pdb_id = "1MBO"
pdb_list = PDBList()
pdb_path = pdb_list.retrieve_pdb_file(pdb_id, pdir='.', file_format="pdb")
pdb_file = pdb_path.replace('.ent', '.pdb')
os.rename(pdb_path, pdb_file)

parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdb_file)

# obliczanie kątów phi i psi
phi_psi_data = calculate_phi_psi_angles(structure, '/usr/bin/dssp')

# wykres Ramachandrana
plot_ramachandran(phi_psi_data)
