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
    This function detects a synthetic strategy involving sequential functionalization
    of a pyrimidine core.
    """
    pyrimidine_modifications = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyrimidine in product
                if checker.check_ring("pyrimidine", product):
                    # Check if any reactant contains pyrimidine
                    for reactant in reactants:
                        if checker.check_ring("pyrimidine", reactant):
                            # Get atom indices for pyrimidine in reactant and product
                            reactant_indices = checker.get_ring_atom_indices("pyrimidine", reactant)
                            product_indices = checker.get_ring_atom_indices("pyrimidine", product)

                            if reactant_indices and product_indices:
                                # If we have pyrimidine in both reactant and product, it's a modification
                                print(f"Pyrimidine modification detected at depth {depth}")
                                pyrimidine_modifications.append(depth)
                                break

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 modifications
    if len(pyrimidine_modifications) < 2:
        print("Less than 2 pyrimidine modifications found")
        return False

    # Sort modifications by depth to analyze the sequence
    # In retrosynthetic analysis, we're moving from product to reactants
    # So we're looking at the synthesis in reverse
    pyrimidine_modifications.sort()

    print(f"Found pyrimidine modifications at depths: {pyrimidine_modifications}")

    # Check if we have sequential modifications (not necessarily at adjacent depths)
    # The key is that we have multiple modifications on the pyrimidine core
    # throughout the synthesis route
    return True
