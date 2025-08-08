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
    This function detects if the synthesis involves formation of a guanidine or urea group.
    """
    has_formation = False

    def dfs_traverse(node):
        nonlocal has_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check for urea formation reactions directly
            if (
                checker.check_reaction("Urea synthesis via isocyanate and primary amine", rsmi)
                or checker.check_reaction("Urea synthesis via isocyanate and secondary amine", rsmi)
                or checker.check_reaction("Urea synthesis via isocyanate and diazo", rsmi)
                or checker.check_reaction("Urea synthesis via isocyanate and sulfonamide", rsmi)
                or checker.check_reaction("urea", rsmi)
            ):
                print(f"Detected urea formation reaction: {rsmi}")
                has_formation = True
                return

            # Check for urea in product but not in reactants
            product_has_urea = checker.check_fg("Urea", product_part)
            reactants_have_urea = any(checker.check_fg("Urea", r) for r in reactants if r)

            if product_has_urea and not reactants_have_urea:
                print(f"Detected urea formation: {rsmi}")
                has_formation = True
                return

            # Check for thiourea formation
            product_has_thiourea = checker.check_fg("Thiourea", product_part)
            reactants_have_thiourea = any(checker.check_fg("Thiourea", r) for r in reactants if r)

            if product_has_thiourea and not reactants_have_thiourea:
                print(f"Detected thiourea formation: {rsmi}")
                has_formation = True
                return

            # Check for isocyanate reacting with amine
            has_isocyanate = any(checker.check_fg("Isocyanate", r) for r in reactants if r)
            has_amine = any(
                checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                for r in reactants
                if r
            )

            if has_isocyanate and has_amine and (product_has_urea or product_has_thiourea):
                print(f"Detected isocyanate-amine reaction to form urea/thiourea: {rsmi}")
                has_formation = True
                return

        # Traverse children
        for child in node.get("children", []):
            if not has_formation:  # Stop traversal if we already found what we're looking for
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Guanidine/urea formation strategy: {has_formation}")
    return has_formation
