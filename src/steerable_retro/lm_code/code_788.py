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
    This function detects a strategy where Boc protection of a nitrogen
    (typically in a piperidine) occurs as the final step of the synthesis.
    """
    found_boc_protection = False
    min_boc_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_boc_protection, min_boc_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants_smiles = parts[0].split(".")
            product_smiles = parts[2]

            print(f"Checking reaction at depth {depth}")
            print(f"Reactants: {reactants_smiles}")
            print(f"Product: {product_smiles}")

            # Check if this is a Boc protection reaction
            is_boc_protection = (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            )

            print(f"Direct reaction check result: {is_boc_protection}")

            # If reaction check fails, try to determine by functional group changes
            if not is_boc_protection:
                # Check if product contains Boc group
                has_boc_in_product = checker.check_fg("Boc", product_smiles)
                print(f"Product has Boc: {has_boc_in_product}")

                if has_boc_in_product:
                    # Check if any reactant has an unprotected nitrogen
                    for reactant in reactants_smiles:
                        has_boc_in_reactant = checker.check_fg("Boc", reactant)
                        print(f"Reactant has Boc: {has_boc_in_reactant}")

                        # Check for various types of amines
                        has_amine = (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        )

                        # Check for piperidine specifically
                        has_piperidine = checker.check_ring("piperidine", reactant)

                        print(f"Reactant has amine: {has_amine}")
                        print(f"Reactant has piperidine: {has_piperidine}")

                        if (has_amine or has_piperidine) and not has_boc_in_reactant:
                            print(f"Found Boc protection at depth {depth} (FG analysis)")
                            is_boc_protection = True
                            break

            if is_boc_protection:
                print(f"Found Boc protection at depth {depth}")
                found_boc_protection = True
                min_boc_depth = min(min_boc_depth, depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Minimum Boc protection depth: {min_boc_depth}")
    # Consider it late-stage if it occurs at depth 0, 1, 2, or 3
    return found_boc_protection and min_boc_depth <= 3
