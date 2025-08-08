#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects if benzyl ether groups are preserved throughout the synthesis.
    """
    # Track benzyl ethers and their depths
    benzyl_ethers_by_depth = {}

    def is_benzyl_ether(mol_smiles):
        """Check if molecule contains benzyl ether groups"""
        # Check for both ether functional group and benzene ring
        if not checker.check_fg("Ether", mol_smiles):
            print(f"No ether found in {mol_smiles}")
            return False

        if not checker.check_ring("benzene", mol_smiles):
            print(f"No benzene found in {mol_smiles}")
            return False

        # Get the molecule
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            print(f"Could not parse molecule: {mol_smiles}")
            return False

        # Get ether oxygen atoms
        ether_matches = checker.get_fg_atom_indices("Ether", mol_smiles)
        if not ether_matches:
            print(f"No ether atom indices found in {mol_smiles}")
            return False

        # Get benzene ring atoms
        benzene_matches = checker.get_ring_atom_indices("benzene", mol_smiles)
        if not benzene_matches:
            print(f"No benzene atom indices found in {mol_smiles}")
            return False

        print(f"Ether matches: {ether_matches}")
        print(f"Benzene matches: {benzene_matches}")

        # Check for benzyl ether pattern: benzene-CH2-O-R
        for ether_match in ether_matches:
            for benzene_match in benzene_matches:
                # Flatten the matches if they're nested
                if isinstance(ether_match[0], tuple):
                    ether_atoms = [atom for sublist in ether_match for atom in sublist]
                else:
                    ether_atoms = ether_match

                if isinstance(benzene_match[0], tuple):
                    benzene_atoms = [atom for sublist in benzene_match for atom in sublist]
                else:
                    benzene_atoms = benzene_match

                # Find oxygen atoms in the ether match
                for atom_idx in ether_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == "O":
                        # Check neighbors of oxygen
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetSymbol() == "C":
                                # Check if this carbon is connected to benzene
                                for benzene_atom_idx in benzene_atoms:
                                    benzene_atom = mol.GetAtomWithIdx(benzene_atom_idx)
                                    for benzene_neighbor in benzene_atom.GetNeighbors():
                                        if benzene_neighbor.GetIdx() == neighbor.GetIdx():
                                            # Found a carbon connected to both oxygen and benzene
                                            print(f"Found benzyl ether pattern in {mol_smiles}")
                                            return True

        print(f"No benzyl ether pattern found in {mol_smiles}")
        return False

    def dfs_traverse(node, depth=0):
        """Traverse the synthesis route and track benzyl ethers by depth"""
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for benzyl ether in this molecule
            if is_benzyl_ether(mol_smiles):
                print(f"Found benzyl ether at depth {depth}: {mol_smiles}")
                if depth not in benzyl_ethers_by_depth:
                    benzyl_ethers_by_depth[depth] = []
                benzyl_ethers_by_depth[depth].append(mol_smiles)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if benzyl ethers are present in at least 3 different depths
    result = len(benzyl_ethers_by_depth) >= 3

    print(f"Benzyl ether preservation detected: {result}")
    print(f"Found benzyl ethers at {len(benzyl_ethers_by_depth)} different depths")
    for depth, smiles_list in benzyl_ethers_by_depth.items():
        print(f"Depth {depth}: {len(smiles_list)} benzyl ethers")

    return result
