#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a synthetic strategy where a nitro group is reduced to an amine
    in the final step of the synthesis.
    """
    final_reaction_has_nitro_reduction = False
    target_molecule = route["smiles"]  # The target molecule is at the root of the route

    print(f"Target molecule: {target_molecule}")

    def dfs_traverse(node, depth=0):
        nonlocal final_reaction_has_nitro_reduction

        print(f"Examining node at depth {depth}, type: {node['type']}")

        # Process reaction nodes
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Reaction at depth {depth}: {rsmi}")

            # Check if this is the first reaction in retrosynthesis (last in forward synthesis)
            # This corresponds to depth 1 in the retrosynthetic tree
            if depth == 1:
                print(f"Checking if this is the final reaction (depth={depth})")

                # Primary check: Use the specific reaction checker
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(
                        "Confirmed nitro reduction in final step using reaction checker"
                    )
                    final_reaction_has_nitro_reduction = True
                else:
                    # Secondary check: In retrosynthesis, the reactant has nitro and product has amine
                    # Note: In retrosynthesis, reactants are precursors and product is target
                    has_nitro_reactant = any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    )
                    has_amine_product = checker.check_fg("Primary amine", product)

                    print(
                        f"Nitro in reactants: {has_nitro_reactant}, Amine in product: {has_amine_product}"
                    )

                    if has_nitro_reactant and has_amine_product:
                        print("Detected nitro to amine conversion in final step")
                        final_reaction_has_nitro_reduction = True

        # Continue traversing the synthesis route
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {final_reaction_has_nitro_reduction}")
    return final_reaction_has_nitro_reduction
