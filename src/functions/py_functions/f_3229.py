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
    Detects reduction of carboxylic acid to aldehyde in the synthesis.
    """
    reduction_found = False

    def dfs_traverse(node):
        nonlocal reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # In forward reaction: reactants -> product
            # In retrosynthesis: product <- reactants
            reactants_forward = rsmi.split(">")[0].split(".")
            product_forward = rsmi.split(">")[-1]

            # Check if this is a known oxidation reaction (which in retrosynthesis would be a reduction)
            if checker.check_reaction(
                "Oxidation of aldehydes to carboxylic acids", rsmi
            ):
                print(
                    f"Found oxidation of aldehyde to carboxylic acid reaction: {rsmi}"
                )
                reduction_found = True
                return

            # Check for specific reduction reactions
            # In retrosynthesis, we're looking at the reverse of these reactions
            reduction_reactions = [
                "Reduction of carboxylic acid to aldehyde",
                "Oxidation of aldehyde to carboxylic acid",
            ]

            for rxn_name in reduction_reactions:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"Found {rxn_name} reaction: {rsmi}")
                    reduction_found = True
                    return

            # If not a known reaction, check for the functional group transformation
            # In forward direction: check if carboxylic acid in reactants and aldehyde in product
            acid_in_reactants = False
            aldehyde_in_product = False

            for reactant in reactants_forward:
                if checker.check_fg("Carboxylic acid", reactant):
                    acid_in_reactants = True
                    acid_reactant = reactant
                    print(f"Found carboxylic acid in reactant: {reactant}")
                    break

            if acid_in_reactants and checker.check_fg("Aldehyde", product_forward):
                aldehyde_in_product = True
                print(f"Found aldehyde in product: {product_forward}")

            if acid_in_reactants and aldehyde_in_product:
                # This is a potential carboxylic acid reduction to aldehyde
                print(
                    "Detected potential carboxylic acid reduction to aldehyde based on functional groups"
                )
                reduction_found = True
                return

            # Check the reverse direction (for retrosynthesis)
            # In retrosynthesis: check if aldehyde in reactants and carboxylic acid in product
            aldehyde_in_reactants = False

            for reactant in reactants_forward:
                if checker.check_fg("Aldehyde", reactant):
                    aldehyde_in_reactants = True
                    print(f"Found aldehyde in reactant: {reactant}")
                    break

            if aldehyde_in_reactants and checker.check_fg(
                "Carboxylic acid", product_forward
            ):
                print(f"Found carboxylic acid in product: {product_forward}")
                print(
                    "Detected oxidation of aldehyde to carboxylic acid (reverse of target transformation)"
                )
                reduction_found = True
                return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if reduction_found:
        print("Detected carboxylic acid reduction to aldehyde strategy")
        return True
    else:
        print("No carboxylic acid reduction to aldehyde detected")
        return False


# --------------------------------------------------
