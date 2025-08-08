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
    Detects if the synthesis includes an amine to nitro transformation.
    In retrosynthetic direction, this means detecting nitro reduction to amine.
    """
    transformation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal transformation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # In retrosynthesis, we're looking for nitro reduction to amine
            # So check for nitro in reactants and amine in product

            # First check if this is a nitro reduction reaction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Found nitro reduction reaction at depth {depth}: {rsmi}")
                transformation_found = True
                return

            # If specific reaction check fails, check for the functional groups
            nitro_in_reactants = False
            for r_smi in reactants_smiles:
                if checker.check_fg("Nitro group", r_smi):
                    nitro_in_reactants = True
                    break

            amine_in_product = checker.check_fg("Primary amine", product_smiles)

            if nitro_in_reactants and amine_in_product:
                # Additional verification could be done here with atom mapping
                # to ensure the nitro group is converted to the amine
                print(f"Found potential nitro to amine transformation at depth {depth}")
                transformation_found = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return transformation_found
