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
    Detects if the synthesis uses Boc protection of an amine group.
    """
    found_boc_protection = False

    def dfs_traverse(node):
        nonlocal found_boc_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a Boc protection or deprotection reaction
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                or checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
            ):
                print(f"Found Boc protection/deprotection reaction: {rsmi}")
                found_boc_protection = True

            # Alternative check: look for carbamic ester (Boc group) in reactants or products
            else:
                # Check for carbamic ester (Boc group)
                carbamic_ester_in_reactants = any(
                    checker.check_fg("Carbamic ester", reactant) for reactant in reactants
                )
                carbamic_ester_in_product = checker.check_fg("Carbamic ester", product)

                # Check for amine groups
                amine_in_reactants = any(
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    for reactant in reactants
                )

                amine_in_product = checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                )

                # If carbamic ester appears in product but not in reactants, and there's an amine in reactants
                if (
                    carbamic_ester_in_product
                    and not carbamic_ester_in_reactants
                    and amine_in_reactants
                ):
                    print(f"Found Boc protection (FG analysis): {rsmi}")
                    found_boc_protection = True

                # If carbamic ester disappears from reactants to product, and there's an amine in product
                elif (
                    carbamic_ester_in_reactants
                    and not carbamic_ester_in_product
                    and amine_in_product
                ):
                    print(f"Found Boc deprotection (FG analysis): {rsmi}")
                    found_boc_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_boc_protection
