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
    Detects a strategy involving protection and subsequent deprotection
    of an indazole N-H, particularly with a tetrahydropyran group.
    """
    protection_found = False
    deprotection_found = False
    protection_depth = float("inf")
    deprotection_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found, protection_depth, deprotection_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check for indazole N-protection (N-H to N-THP)
                reactants_have_indazole = any(checker.check_ring("indazole", r) for r in reactants)
                product_has_indazole = checker.check_ring("indazole", product)
                reactants_have_thp = any(
                    checker.check_ring("tetrahydropyran", r) for r in reactants
                )
                product_has_thp = checker.check_ring("tetrahydropyran", product)

                # Protection: Indazole in reactants, indazole+THP in product
                if (
                    reactants_have_indazole
                    and product_has_indazole
                    and not reactants_have_thp
                    and product_has_thp
                ):
                    # Check if this is likely an N-protection (not perfect but reasonable heuristic)
                    protection_found = True
                    protection_depth = min(protection_depth, depth)
                    print(f"Indazole N-protection detected at depth {depth}")
                    print(f"  Reaction: {rsmi}")

                # Deprotection: Indazole+THP in reactants, indazole (no THP) in product
                if (
                    reactants_have_indazole
                    and product_has_indazole
                    and reactants_have_thp
                    and not product_has_thp
                ):
                    deprotection_found = True
                    deprotection_depth = min(deprotection_depth, depth)
                    print(f"Indazole N-deprotection detected at depth {depth}")
                    print(f"  Reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Protection found: {protection_found} at depth {protection_depth}")
    print(f"Deprotection found: {deprotection_found} at depth {deprotection_depth}")

    # Ensure protection happens before deprotection (higher depth number)
    # In retrosynthetic analysis, protection should be at higher depth than deprotection
    return protection_found and deprotection_found and protection_depth > deprotection_depth
