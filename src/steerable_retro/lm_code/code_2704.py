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
    This function detects a synthetic strategy involving a persistent benzyloxy group
    throughout the synthesis.
    """
    # Track molecules that contain benzyloxy groups
    molecules_with_benzyloxy = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check if the molecule contains a benzyloxy group
            if checker.check_fg("Ether", mol_smiles):
                # More specific check for benzyloxy pattern
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    benzyloxy_pattern = Chem.MolFromSmarts("c[CH2][O]")
                    if mol.HasSubstructMatch(benzyloxy_pattern):
                        print(f"Found benzyloxy group in molecule: {mol_smiles}")
                        molecules_with_benzyloxy.append((mol_smiles, depth))

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # For reaction nodes, we don't need to do anything specific
            pass

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found benzyloxy groups and if they persist throughout the synthesis
    if not molecules_with_benzyloxy:
        print("No benzyloxy groups found in the synthesis")
        return False

    # Sort by depth to check persistence from early to late stage
    molecules_with_benzyloxy.sort(key=lambda x: x[1], reverse=True)

    # Check if benzyloxy is present in at least 3 consecutive steps
    if len(molecules_with_benzyloxy) >= 3:
        print(
            f"Found persistent benzyloxy group throughout synthesis in {len(molecules_with_benzyloxy)} molecules"
        )
        return True
    else:
        print(
            f"Benzyloxy group not persistent enough, only found in {len(molecules_with_benzyloxy)} molecules"
        )
        return False
