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
    This function detects a synthetic strategy involving epoxide ring-opening
    by a nucleophilic amine.
    """
    has_epoxide_opening = False

    def dfs_traverse(node, depth=0):
        nonlocal has_epoxide_opening

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Direct check for epoxide ring opening with amine reaction
                if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                    print(f"Detected epoxide ring-opening reaction at depth {depth}")
                    has_epoxide_opening = True
                else:
                    # Fallback check using functional groups
                    # Check for epoxide in reactants
                    has_epoxide_reactant = False
                    for reactant in reactants_smiles:
                        if checker.check_ring("oxirane", reactant):
                            has_epoxide_reactant = True
                            print(f"Found epoxide in reactant: {reactant}")
                            break

                    # Check for amine in reactants
                    has_amine_reactant = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            has_amine_reactant = True
                            print(f"Found amine in reactant: {reactant}")
                            break

                    # Check for alcohol in product (result of epoxide opening)
                    has_alcohol_product = checker.check_fg(
                        "Primary alcohol", product_smiles
                    ) or checker.check_fg("Secondary alcohol", product_smiles)
                    if has_alcohol_product:
                        print(f"Found alcohol in product: {product_smiles}")

                    # Check for amine in product (should be retained from reactant)
                    has_amine_product = checker.check_fg(
                        "Primary amine", product_smiles
                    ) or checker.check_fg("Secondary amine", product_smiles)
                    if has_amine_product:
                        print(f"Found amine in product: {product_smiles}")

                    if (
                        has_epoxide_reactant
                        and has_amine_reactant
                        and has_alcohol_product
                        and has_amine_product
                    ):
                        print(f"Detected epoxide ring-opening by pattern matching at depth {depth}")
                        has_epoxide_opening = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_epoxide_opening
