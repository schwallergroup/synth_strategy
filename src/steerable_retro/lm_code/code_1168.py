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
    Detects if the synthesis route includes sequential nitrogen functional group interconversions.
    """
    # Track nitrogen functional group transformations
    n_transformations = []

    def dfs(node, depth=0):
        nonlocal n_transformations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for various nitrogen functional group interconversions
            nitrogen_reactions = [
                "Reduction of nitro groups to amines",
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Boc amine protection",
                "Boc amine deprotection",
            ]

            for reaction_type in nitrogen_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    # Verify nitrogen functional group change
                    n_fg_before = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        or checker.check_fg("Nitro group", r)
                        or checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        for r in reactants
                    )

                    n_fg_after = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                    )

                    if n_fg_before and n_fg_after:
                        n_transformations.append((reaction_type, depth))
                        print(f"Found nitrogen transformation: {reaction_type} at depth {depth}")
                    break

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    # Check if we have at least two sequential nitrogen transformations
    # Sort by depth to check if they're sequential
    n_transformations.sort(key=lambda x: x[1])

    # Need at least two transformations to be sequential (adjacent depths)
    return len(n_transformations) >= 2 and any(
        abs(n_transformations[i][1] - n_transformations[i + 1][1]) == 1
        for i in range(len(n_transformations) - 1)
    )
