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
    Detects if the final step in the synthesis involves reduction of a nitro group to an amine
    """
    # Track if we found nitro reduction at depth 1 (final synthetic step)
    found_nitro_reduction_at_final_step = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction_at_final_step

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Final step in synthesis is at depth 1 (first reaction from target)
            if depth == 1:
                print(f"Analyzing final step reaction: {rsmi}")

                # Primary check: Use the reaction checker
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print("Found nitro reduction at final step using reaction checker")
                    found_nitro_reduction_at_final_step = True
                    return

                # Secondary check: Look for specific nitro to amine conversion
                # We need to ensure the nitro group is actually being reduced to an amine
                # in the same molecule structure

                # Check product for amine groups
                has_amine = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                if has_amine:
                    print("Product contains amine group")

                    # Check reactants for nitro groups
                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant):
                            print(f"Reactant contains nitro group: {reactant}")

                            # Create molecule objects for more detailed analysis
                            try:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                product_mol = Chem.MolFromSmiles(product)

                                if reactant_mol and product_mol:
                                    # Use atom mapping to verify the transformation
                                    # In a nitro reduction, the core structure should remain similar
                                    # while the nitro group is converted to an amine

                                    # Check if this is a reduction reaction (common reagents)
                                    reduction_reagents = [
                                        "[H]",
                                        "H2",
                                        "Pd",
                                        "Pt",
                                        "Ni",
                                        "Fe",
                                        "Zn",
                                        "Sn",
                                    ]
                                    has_reduction_reagent = any(
                                        re.search(reagent, rsmi)
                                        for reagent in reduction_reagents
                                    )

                                    if has_reduction_reagent:
                                        print("Reaction contains reduction reagents")
                                        found_nitro_reduction_at_final_step = True
                                        return
                            except:
                                print("Error in detailed molecule analysis")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_nitro_reduction_at_final_step}")

    return found_nitro_reduction_at_final_step
