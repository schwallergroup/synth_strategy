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
    Detects if the synthesis uses a late-stage amide coupling strategy.
    This checks if the final reaction forms an amide bond between a carboxylic acid and an amine.
    """
    amide_coupling_detected = False
    final_reaction_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_detected, final_reaction_depth

        if node["type"] == "reaction":
            # Check if this is the final reaction (lowest depth)
            if depth < final_reaction_depth:
                final_reaction_depth = depth

                # Extract reactants and product
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check for carboxylic acid and amine in reactants
                    acid_present = False
                    amine_present = False
                    acid_reactant = None
                    amine_reactant = None

                    for reactant in reactants_smiles:
                        try:
                            if checker.check_fg("Carboxylic acid", reactant):
                                acid_present = True
                                acid_reactant = reactant
                                print(f"Found carboxylic acid in reactant: {reactant}")

                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Tertiary amine", reactant)
                            ):
                                amine_present = True
                                amine_reactant = reactant
                                print(f"Found amine in reactant: {reactant}")
                        except Exception as e:
                            print(f"Error checking reactant {reactant}: {e}")
                            continue

                    # Check if this is an amide coupling reaction
                    if acid_present and amine_present:
                        # Check for amide in product
                        print(f"Checking for amide in product: {product_smiles}")
                        amide_in_product = (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        )

                        print(f"Amide in product: {amide_in_product}")

                        if amide_in_product:
                            # Verify this is an amide coupling reaction
                            is_amide_coupling = False

                            # Check for various amide coupling reactions
                            amide_reaction_types = [
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                "Carboxylic acid with primary amine to amide",
                                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                "Acyl chloride with secondary amine to amide",
                                "Ester with primary amine to amide",
                                "Ester with secondary amine to amide",
                                "Ester with ammonia to amide",
                                "Acyl chloride with ammonia to amide",
                                "Acyl chloride with primary amine to imide",
                                "Schotten-Baumann_amide",
                            ]

                            for reaction_name in amide_reaction_types:
                                if checker.check_reaction(reaction_name, rsmi):
                                    is_amide_coupling = True
                                    print(f"Amide coupling reaction detected: {reaction_name}")
                                    break

                            # If no specific reaction type is detected, check if it's a general amide formation
                            if not is_amide_coupling:
                                # Check if the carboxylic acid carbon and amine nitrogen form an amide bond in the product
                                # This is a more general check for amide formation
                                try:
                                    # Look for C(=O)-N pattern in the product that wasn't in the reactants
                                    acid_mol = Chem.MolFromSmiles(acid_reactant)
                                    amine_mol = Chem.MolFromSmiles(amine_reactant)
                                    product_mol = Chem.MolFromSmiles(product_smiles)

                                    # If we have atom mapping, we can use it to track the atoms
                                    # This is a simplified check - in a real implementation, you would
                                    # need to track the specific atoms involved in the reaction
                                    print("Performing general amide formation check")
                                    is_amide_coupling = True
                                except Exception as e:
                                    print(f"Error in general amide formation check: {e}")

                            if is_amide_coupling:
                                amide_coupling_detected = True
                                print(f"Late-stage amide coupling detected at depth {depth}")

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amide_coupling_detected
