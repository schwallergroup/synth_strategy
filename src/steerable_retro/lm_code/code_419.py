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
    Detects if the synthesis route involves a late-stage esterification
    (carboxylic acid to ester transformation in the final step).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        print(f"Examining node of type {node['type']} at depth {depth}")

        if node["type"] == "reaction" and depth <= 1:  # Check final or near-final reactions
            print(f"Examining potential late-stage reaction at depth {depth}")

            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # In retrosynthesis, the product is what we're breaking down
                # So in the forward direction, the product is the target molecule
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is an esterification reaction using multiple reaction types
                is_esterification = (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                )

                if is_esterification:
                    print(f"Found esterification reaction: {rsmi}")
                    result = True
                    return

                # For retrosynthesis, we need to check if:
                # 1. The product (target molecule) contains an ester
                # 2. The reactants contain a carboxylic acid and/or alcohol

                # Check if product has an ester group
                ester_in_product = checker.check_fg("Ester", product_smiles)
                if ester_in_product:
                    print(f"Found ester in product: {product_smiles}")

                    # Check reactants for carboxylic acid and alcohol
                    carboxylic_acid_found = False
                    alcohol_found = False
                    acyl_halide_found = False

                    for reactant in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant):
                            carboxylic_acid_found = True
                            print(f"Found carboxylic acid in reactant: {reactant}")

                        # Check for any type of alcohol
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                            or checker.check_fg("Phenol", reactant)
                        ):
                            alcohol_found = True
                            print(f"Found alcohol in reactant: {reactant}")

                        if checker.check_fg("Acyl halide", reactant):
                            acyl_halide_found = True
                            print(f"Found acyl halide in reactant: {reactant}")

                    # If we have carboxylic acid and alcohol in reactants, and ester in product,
                    # it's likely an esterification
                    if (carboxylic_acid_found and alcohol_found) or (
                        acyl_halide_found and alcohol_found
                    ):
                        print(
                            "Found late-stage esterification based on functional group transformation"
                        )
                        result = True
                        return

                # Also check for the reverse direction in retrosynthesis
                # In retrosynthesis, we might be breaking down an ester into acid + alcohol
                ester_in_reactants = any(
                    checker.check_fg("Ester", reactant) for reactant in reactants_smiles
                )
                if ester_in_reactants:
                    print(f"Found ester in reactants")

                    # Check if product has carboxylic acid
                    if checker.check_fg("Carboxylic acid", product_smiles):
                        print(f"Found carboxylic acid in product: {product_smiles}")
                        print("Found late-stage esterification (hydrolysis in retrosynthesis)")
                        result = True
                        return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {result}")
    return result
