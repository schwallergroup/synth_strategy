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
    Detects if there's retention of aromatic bromide through multiple steps.
    """
    # Track molecules with aromatic bromides at different depths
    bromide_molecules = []
    bromide_retention = False

    def dfs(node, depth=0):
        nonlocal bromide_retention

        if node.get("type") == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check if molecule has aromatic bromide
            if checker.check_fg("Aromatic halide", mol_smiles) and "Br" in mol_smiles:
                bromide_molecules.append((depth, mol_smiles))
                print(f"Aromatic bromide found at depth {depth}: {mol_smiles}")

                # Check if we have bromides at different depths
                depths = set(d for d, _ in bromide_molecules)
                if len(depths) >= 2:
                    bromide_retention = True
                    print(f"Aromatic bromide retention detected across depths: {sorted(depths)}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS from the root
    dfs(route)

    return bromide_retention
