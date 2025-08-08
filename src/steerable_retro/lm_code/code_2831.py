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
    This function detects if the synthetic route involves dichlorinated aromatic compounds.
    """
    dichlorinated_aromatic_found = False

    def dfs_traverse(node):
        nonlocal dichlorinated_aromatic_found

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # Check if molecule has aromatic halide functional group
                if checker.check_fg("Aromatic halide", mol_smiles):
                    # Count chlorines attached to aromatic rings
                    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]

                    # Count chlorines attached to aromatic atoms
                    chlorine_count = 0
                    for atom in aromatic_atoms:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetSymbol() == "Cl":
                                chlorine_count += 1

                    # Check if there are exactly 2 chlorines attached to aromatic rings
                    if chlorine_count == 2:
                        dichlorinated_aromatic_found = True
                        print(f"Found dichlorinated aromatic: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Dichlorinated aromatic strategy detected: {dichlorinated_aromatic_found}")
    return dichlorinated_aromatic_found
