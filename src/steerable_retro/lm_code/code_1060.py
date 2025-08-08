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
    This function detects a synthetic strategy where one heterocycle (indole)
    is maintained while another heterocycle (thiazole) is constructed.
    """
    # Track if we found the pattern
    found_indole_in_starting_materials = False
    found_thiazole_in_starting_materials = False
    found_indole_thiazole_in_product = False
    thiazole_constructed = False
    final_product_smiles = None

    # First pass to identify the final product (leaf molecule node)
    def find_final_product(node, depth=0):
        nonlocal final_product_smiles

        if node["type"] == "mol" and not node.get("in_stock", False) and not node.get("children"):
            # This is a leaf molecule node that's not a starting material - potential final product
            if final_product_smiles is None or depth < find_final_product.min_depth:
                final_product_smiles = node["smiles"]
                find_final_product.min_depth = depth

        # Traverse children
        for child in node.get("children", []):
            find_final_product(child, depth + 1)

    find_final_product.min_depth = float("inf")
    find_final_product(route)

    # If no final product found, try to find the molecule at the lowest depth
    if final_product_smiles is None:

        def find_lowest_depth_molecule(node, depth=0):
            nonlocal final_product_smiles

            if node["type"] == "mol":
                if final_product_smiles is None or depth < find_lowest_depth_molecule.min_depth:
                    final_product_smiles = node["smiles"]
                    find_lowest_depth_molecule.min_depth = depth

            # Traverse children
            for child in node.get("children", []):
                find_lowest_depth_molecule(child, depth + 1)

        find_lowest_depth_molecule.min_depth = float("inf")
        find_lowest_depth_molecule(route)

    print(f"Identified final product: {final_product_smiles}")

    # Second pass to check for the pattern
    def dfs_traverse(node, depth=0):
        nonlocal found_indole_in_starting_materials, found_thiazole_in_starting_materials
        nonlocal found_indole_thiazole_in_product, thiazole_constructed

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_indole = checker.check_ring("indole", mol_smiles)
            has_thiazole = checker.check_ring("thiazole", mol_smiles)

            # Check if this is a starting material
            if node.get("in_stock", False):
                # Check if starting material contains indole
                if has_indole:
                    found_indole_in_starting_materials = True
                    print(f"Found indole in starting material: {mol_smiles}")

                # Check if starting material contains thiazole
                if has_thiazole:
                    found_thiazole_in_starting_materials = True
                    print(f"Found thiazole in starting material: {mol_smiles}")

            # Check if this is the final product
            if mol_smiles == final_product_smiles:
                if has_indole and has_thiazole:
                    found_indole_thiazole_in_product = True
                    print(f"Found both indole and thiazole in final product: {mol_smiles}")
                    print(f"Indole present: {has_indole}, Thiazole present: {has_thiazole}")

        elif node["type"] == "reaction":
            # Check if this reaction constructs a thiazole ring
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has thiazole but reactants don't
                product_has_thiazole = checker.check_ring("thiazole", product)
                reactants_have_thiazole = any(
                    checker.check_ring("thiazole", reactant) for reactant in reactants
                )

                if product_has_thiazole and not reactants_have_thiazole:
                    thiazole_constructed = True
                    print(f"Thiazole ring constructed in reaction: {rsmi}")

                # Also check if the reaction is a known thiazole-forming reaction
                if checker.check_reaction("thiazole", rsmi) or checker.check_reaction(
                    "{thiazole}", rsmi
                ):
                    thiazole_constructed = True
                    print(f"Detected thiazole-forming reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # The pattern is found if:
    # 1. Indole is present in starting materials
    # 2. Thiazole is NOT present in starting materials (it's constructed)
    # 3. Both indole and thiazole are present in the final product
    # 4. A thiazole ring is constructed during synthesis
    result = (
        found_indole_in_starting_materials
        and not found_thiazole_in_starting_materials
        and found_indole_thiazole_in_product
        and thiazole_constructed
    )

    print(f"Indole in starting materials: {found_indole_in_starting_materials}")
    print(f"Thiazole in starting materials: {found_thiazole_in_starting_materials}")
    print(f"Both indole and thiazole in final product: {found_indole_thiazole_in_product}")
    print(f"Thiazole constructed during synthesis: {thiazole_constructed}")
    print(f"Final result: {result}")

    return result
