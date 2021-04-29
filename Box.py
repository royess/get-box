from Bio.PDB import *
import numpy as np

def nearest_dist(residue : Residue.Residue, ligand : Residue.Residue):
    dists = [[(a.get_vector() - aa.get_vector()).norm() for a in residue] for aa in ligand]
    return np.min(dists)


def get_residues_in_pocket(residues, ligand, cutoff=3):
    ligand_id = ligand.get_id()
    nearest_dists = [nearest_dist(r, ligand) for r in residues if r.get_resname() not in ['HOH'] if r.get_id()!=ligand_id]
    idxs = [i for i in range(len(nearest_dists)) if nearest_dists[i] < cutoff]
    return [list(residues)[i] for i in idxs]


def get_box(residues, add=4.0):
    coords = []
    for r in residues:
        coords += [list(a.get_vector()) for a in r]
    coords = np.array(coords)
    vecu = coords.max(axis=0)
    vecl = coords.min(axis=0)
    center = (vecu + vecl) / 2.0
    sizes = vecu - vecl + add
    return center, sizes


def get_box_config(file_path, residues, ligand, cutoff=3, add_size=4, verbose=True):
    residues_in_pocket = get_residues_in_pocket(residues, ligand, cutoff)
    if(verbose):
        print("Residues in pocket:")
        for r in residues_in_pocket:
            print(r)
    
    center, sizes = get_box(residues_in_pocket)
    text = f'''\
    center_x = {center[0]}
    center_y = {center[1]}
    center_z = {center[2]}
    size_x = {sizes[0] + add_size}
    size_y = {sizes[1] + add_size}
    size_z = {sizes[2] + add_size}'''
    
    if(verbose):
        print("Box:")
        print(text)
    with open(file_path, 'w+') as f:
        f.write(text)
    
