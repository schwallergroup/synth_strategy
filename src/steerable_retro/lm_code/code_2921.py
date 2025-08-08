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
    Detects if the synthesis involves cleavage of an aryl ether.
    """
    ether_cleavage_found = False

    def dfs_traverse(node, depth=0):
        nonlocal ether_cleavage_found

        if ether_cleavage_found:
            return  # Early return if we already found what we're looking for

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction SMILES: {rsmi}")

            # In retrosynthetic analysis, the product in rsmi is what we're making
            # and the reactants are what we'd use in the forward direction
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an ether cleavage reaction
            is_ether_cleavage = (
                checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                or checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
            )

            if is_ether_cleavage:
                print(f"Detected ether cleavage reaction by reaction type: {rsmi}")
                ether_cleavage_found = True
                return

            # Even if not detected by reaction type, check for functional group patterns
            # In retrosynthesis, the product contains the ether and reactants contain the alcohol
            product_has_ether = checker.check_fg("Ether", product)

            # Check for various types of alcohols in the reactants
            reactants_have_alcohol = any(
                checker.check_fg("Phenol", r)
                or checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                or checker.check_fg("Aromatic alcohol", r)
                for r in reactants
            )

            if product_has_ether and reactants_have_alcohol:
                print(f"Detected potential ether cleavage by functional group analysis")

                # Additional verification: check for common ether cleavage patterns
                # Look for reactions that might involve ether cleavage but weren't caught by the reaction checkers
                for r in reactants:
                    if (
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Phenol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                    ):
                        print(
                            f"Confirmed ether cleavage: ether in product and alcohol in reactants"
                        )
                        ether_cleavage_found = True
                        return

        # Traverse children (moving backward in synthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ether_cleavage_found
