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
    This function detects a synthetic strategy involving late-stage nitro reduction to form an amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_nitro_reduction

        # For reaction nodes, check if it's a nitro reduction
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Use depth from metadata if available, otherwise use calculated depth
            node_depth = node["metadata"].get("depth", current_depth)

            # Check if this is a late-stage reaction (within first 2 steps)
            if node_depth <= 2:
                print(f"Examining late-stage reaction at depth {node_depth}: {rsmi}")

                # Check if this is a nitro reduction reaction using the reaction checker
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {node_depth}")
                    has_nitro_reduction = True
                else:
                    # Fallback check: verify nitro group in reactants is converted to amine in product
                    nitro_in_reactants = any(
                        checker.check_fg("Nitro group", reactant) for reactant in reactants
                    )
                    amine_in_product = checker.check_fg("Primary amine", product)

                    if nitro_in_reactants and amine_in_product:
                        # Additional check to ensure nitro is actually being reduced
                        # Count nitro groups in reactants and product
                        nitro_count_reactants = sum(
                            1 for reactant in reactants if checker.check_fg("Nitro group", reactant)
                        )
                        nitro_count_product = 1 if checker.check_fg("Nitro group", product) else 0

                        # If nitro groups decreased and amines appeared, it's likely a reduction
                        if nitro_count_reactants > nitro_count_product:
                            print(f"Found nitro to amine conversion at depth {node_depth}")
                            has_nitro_reduction = True

        # Recursively traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_nitro_reduction
