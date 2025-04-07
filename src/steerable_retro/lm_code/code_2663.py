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
    This function detects if the synthesis maintains BOC protection throughout
    the route until the final product.
    """
    # Track BOC protection status for each molecule in the main synthetic pathway
    boc_status = {}
    main_pathway_molecules = set()
    final_product_smiles = route["smiles"]
    main_pathway_molecules.add(final_product_smiles)

    # Track if we've seen a BOC deprotection as the final step
    final_step_is_boc_deprotection = False

    def dfs_traverse(node, depth=0, is_main_pathway=True):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Only track molecules in the main synthetic pathway
            if is_main_pathway:
                main_pathway_molecules.add(mol_smiles)

                # Check if molecule has BOC protection
                has_boc = checker.check_fg("Boc", mol_smiles)
                boc_status[mol_smiles] = has_boc

                # Print for debugging
                print(f"Main pathway molecule at depth {depth}: {mol_smiles}")
                print(f"Has BOC: {has_boc}")

                # If this is a starting material (in_stock=True), we don't require BOC protection
                if node.get("in_stock", False):
                    print(f"Starting material (in_stock=True): {mol_smiles}")
                    # Remove from tracking as it's a starting material
                    if mol_smiles in boc_status:
                        del boc_status[mol_smiles]
            else:
                print(f"Reagent/side molecule at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            # For reaction nodes, check if BOC protection is maintained
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction at depth {depth}: {rsmi}")

                # Check if this is a BOC deprotection reaction
                is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)
                print(f"Is BOC deprotection: {is_boc_deprotection}")

                # If it's a BOC deprotection, we need to check if it's the final step
                if is_boc_deprotection:
                    # Get the product of this reaction
                    product = rsmi.split(">")[-1]

                    # If this product is the final product, mark that we've seen a BOC deprotection as the final step
                    if product == final_product_smiles:
                        print(f"BOC deprotection is the final step at depth {depth}")
                        nonlocal final_step_is_boc_deprotection
                        final_step_is_boc_deprotection = True
                    else:
                        print(f"BOC deprotection occurred too early at depth {depth}")
                        return False

                # For reactions, identify the main reactant for the next step
                # The product of this reaction becomes the main reactant for the next step
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Traverse children with appropriate main pathway flag
                for child in node.get("children", []):
                    if child["type"] == "mol":
                        # If this child is a reactant in the main reaction, it's part of the main pathway
                        child_is_main = (
                            child["smiles"] in reactants
                            and child["smiles"] in main_pathway_molecules
                        )
                        if not dfs_traverse(child, depth + 1, child_is_main):
                            return False
                    else:
                        # For reaction nodes, continue with current main pathway status
                        if not dfs_traverse(child, depth + 1, is_main_pathway):
                            return False

                # Skip the default traversal since we've handled it above
                return True

        # Traverse children
        for child in node.get("children", []):
            if not dfs_traverse(child, depth + 1, is_main_pathway):
                return False

        return True

    # Start traversal
    if not dfs_traverse(route):
        return False

    # Check if final product has BOC protection
    final_has_boc = checker.check_fg("Boc", final_product_smiles)
    print(f"Final product: {final_product_smiles}")
    print(f"Final product has BOC: {final_has_boc}")

    # Check if all intermediates in the main pathway (except possibly the final product) have BOC protection
    all_intermediates_protected = True
    for smiles, has_boc in boc_status.items():
        # Skip the final product
        if smiles == final_product_smiles:
            continue

        if not has_boc:
            print(f"Main pathway intermediate without BOC protection: {smiles}")
            all_intermediates_protected = False
            break

    # The synthesis maintains BOC protection if:
    # 1. All intermediates in the main pathway (except possibly the final product) have BOC protection
    # 2. Either the final product has BOC protection OR the final step is a BOC deprotection
    return all_intermediates_protected
