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
    Detects a strategy involving Boc protection of an indole nitrogen during synthesis.
    """
    # Initialize tracking variables
    has_boc_protection = False
    has_indole = False
    has_boc_protected_indole = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection, has_indole, has_boc_protected_indole

        if node["type"] == "mol":
            if "smiles" in node:
                mol_smiles = node["smiles"]

                # Check for indole
                if checker.check_ring("indole", mol_smiles):
                    has_indole = True
                    print(f"Detected indole at depth {depth}: {mol_smiles}")

                # Check for Boc group and indole together
                if checker.check_fg("Boc", mol_smiles) and checker.check_ring("indole", mol_smiles):
                    has_boc_protected_indole = True
                    has_boc_protection = True
                    print(f"Detected N-Boc protected indole at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for Boc protection reaction
            if checker.check_reaction("Boc amine protection", rsmi):
                has_boc_protection = True
                print(f"Detected Boc protection reaction at depth {depth}: {rsmi}")

            # Also check for Boc deprotection as part of the strategy
            if checker.check_reaction("Boc amine deprotection", rsmi):
                has_boc_protection = True
                print(f"Detected Boc deprotection reaction at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required elements are present
    result = has_boc_protection and has_indole

    print(f"Boc protection detected: {has_boc_protection}")
    print(f"Indole detected: {has_indole}")
    print(f"N-Boc protected indole detected: {has_boc_protected_indole}")
    print(f"Overall Boc protection strategy detected: {result}")

    return result
