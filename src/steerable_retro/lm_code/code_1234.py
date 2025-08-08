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
    Detects if the synthesis involves a nitro reduction in the final step.
    """
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        # Debug information
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Reaction at depth {depth}: {rsmi}")

            # Check if this is a late stage reaction (depth 0 or 1 in retrosynthetic tree)
            if depth <= 1:
                print(f"Checking late stage reaction at depth {depth}: {rsmi}")

                # Check for nitro groups in reactants
                has_nitro = any(checker.check_fg("Nitro group", reactant) for reactant in reactants)
                if has_nitro:
                    print(f"Found nitro group in reactants at depth {depth}")

                # Check for primary amine in product
                has_amine = checker.check_fg("Primary amine", product)
                if has_amine:
                    print(f"Found primary amine in product at depth {depth}")

                # Method 1: Check using the reaction classifier
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction using reaction classifier at depth {depth}")
                    nitro_reduction_found = True
                    return

                # Method 2: Check for nitro group in reactants and primary amine in product
                elif has_nitro and has_amine:
                    print(
                        f"Found nitro group in reactant and primary amine in product at depth {depth}"
                    )
                    nitro_reduction_found = True
                    return

                # Method 3: Check for partial reduction patterns
                elif has_nitro:
                    # If we have a nitro group in reactants, check if it's being reduced
                    if (
                        checker.check_reaction("Reduction", rsmi)
                        or checker.check_reaction("Hydrogenation", rsmi)
                        or checker.check_reaction("Reduction of nitro groups to amines", rsmi)
                    ):
                        print(
                            f"Found potential nitro reduction based on reaction type at depth {depth}"
                        )
                        nitro_reduction_found = True
                        return

        # Continue DFS traversal with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_reduction_found
