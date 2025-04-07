#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects if a chlorophenyl group is maintained throughout the synthesis,
    indicating it's a key structural element that's preserved.
    """
    # Store molecules with their depth in the synthesis route
    mol_data = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            # Only include non-starting materials
            if not node.get("in_stock", False):
                mol_data.append({"smiles": node["smiles"], "depth": depth})

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort by depth (ascending) - final product first, early intermediates last
    mol_data.sort(key=lambda x: x["depth"])

    # Check if chlorophenyl is present in all significant molecules
    complex_mol_threshold = 0  # Include all molecules in the analysis

    complex_mols_with_chlorophenyl = 0
    total_complex_mols = 0

    for mol_info in mol_data:
        smiles = mol_info["smiles"]
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.GetNumAtoms() >= complex_mol_threshold:
                total_complex_mols += 1
                # Check for chlorinated aromatic ring
                if checker.check_fg("Aromatic halide", smiles) and "Cl" in smiles:
                    complex_mols_with_chlorophenyl += 1
                    print(f"Found chlorophenyl in molecule: {smiles}")
        except Exception as e:
            print(f"Error processing molecule {smiles}: {e}")
            continue

    # Consider chlorophenyl maintained if at least 50% of molecules have it
    maintained = (
        total_complex_mols > 0 and complex_mols_with_chlorophenyl >= total_complex_mols * 0.5
    )

    print(
        f"Chlorophenyl maintained: {maintained} ({complex_mols_with_chlorophenyl}/{total_complex_mols})"
    )
    return maintained
