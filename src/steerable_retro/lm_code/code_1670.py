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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects a synthetic strategy where the final step (depth 0)
    involves the reduction of an amide to an amine.
    """
    final_step_is_amide_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_reduction

        # For the root node (final product in retrosynthesis)
        if depth == 0 and node["type"] == "mol":
            print(f"Examining final product: {node['smiles']}")

            # Check immediate children for the final reaction step
            for child in node.get("children", []):
                if (
                    child["type"] == "reaction"
                    and "metadata" in child
                    and "rsmi" in child["metadata"]
                ):
                    rsmi = child["metadata"]["rsmi"]
                    print(f"Examining final step reaction: {rsmi}")

                    # Check if this is an amide reduction reaction using reaction checkers
                    is_amide_reduction = (
                        checker.check_reaction("Reduction of primary amides to amines", rsmi)
                        or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                        or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                        or checker.check_reaction(
                            "Hydrogenolysis of amides/imides/carbamates", rsmi
                        )
                    )

                    if is_amide_reduction:
                        print("Detected amide reduction reaction in final step")
                        final_step_is_amide_reduction = True
                        return

                    # Fallback: Check if reactants have amides and product has amines
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check reactants for amides
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                        ):
                            print(f"Found amide in reactant: {reactant}")

                            # Check if product has amine but no amide
                            product_has_amine = (
                                checker.check_fg("Primary amine", product)
                                or checker.check_fg("Secondary amine", product)
                                or checker.check_fg("Tertiary amine", product)
                            )

                            product_has_amide = (
                                checker.check_fg("Primary amide", product)
                                or checker.check_fg("Secondary amide", product)
                                or checker.check_fg("Tertiary amide", product)
                            )

                            # Verify this is likely an amide reduction
                            if product_has_amine and not product_has_amide:
                                # Try to verify the transformation using atom mapping if available
                                try:
                                    # Check if the reaction involves amide to amine conversion
                                    # by looking at the overall transformation
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    product_mol = Chem.MolFromSmiles(product)

                                    if reactant_mol and product_mol:
                                        # If we have atom mapping, we could check more precisely
                                        # But for now, we'll use the presence of amide and amine
                                        # along with the absence of other major transformations
                                        print(
                                            "Product has amine but no amide - likely an amide reduction"
                                        )
                                        final_step_is_amide_reduction = True
                                        return
                                except Exception as e:
                                    print(f"Error in detailed analysis: {e}")

        # Continue traversal for other nodes
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final result: {final_step_is_amide_reduction}")
    return final_step_is_amide_reduction
