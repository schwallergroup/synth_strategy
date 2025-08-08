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
    This function detects a late-stage nitro reduction strategy where a nitro group
    is reduced to an amine in the final or near-final step of the synthesis.
    """
    has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            try:
                # First check if this is a nitro reduction reaction using the checker
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                # If not found by reaction type, check by functional group transformation
                if not is_nitro_reduction:
                    # Check if reactants have nitro groups and product has amine groups
                    has_nitro_reactant = any(
                        checker.check_fg("Nitro group", reactant) for reactant in reactants_smiles
                    )
                    has_amine_product = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    # Check if nitro groups in reactants but not in product
                    no_nitro_product = not checker.check_fg("Nitro group", product_smiles)

                    # If reactants have nitro, product has amine, and product doesn't have nitro, it's likely a nitro reduction
                    is_nitro_reduction = (
                        has_nitro_reactant and has_amine_product and no_nitro_product
                    )

                if is_nitro_reduction:
                    print(f"Found nitro reduction reaction at depth {depth}")

                    # Check if this is a late-stage reaction (depth 0, 1, or 2)
                    if depth <= 2:
                        print(f"Found late-stage nitro reduction at depth {depth}")
                        has_nitro_reduction = True
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child_idx, child in enumerate(node.get("children", [])):
            # Increment depth when moving from molecule to reaction in retrosynthetic analysis
            child_depth = depth
            if node["type"] == "mol" and child["type"] == "reaction":
                child_depth = depth + 1
            dfs_traverse(child, child_depth)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {has_nitro_reduction}")

    return has_nitro_reduction
