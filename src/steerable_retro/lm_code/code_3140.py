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
    Detects if the route involves late-stage sulfonamide formation (depth 0-2)
    with heterocycle-containing fragments.
    """
    sulfonamide_formed = False
    sulfonamide_depth = -1
    heterocycle_connected_to_sulfonamide = False

    # List of heterocyclic rings to check
    heterocycles = [
        "furan",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "thiophene",
        "benzothiophene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formed, sulfonamide_depth, heterocycle_connected_to_sulfonamide

        if node["type"] == "reaction" and depth <= 2:
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for sulfonamide formation reaction
                is_sulfonamide_reaction = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                # Also check if product contains sulfonamide group
                has_sulfonamide = checker.check_fg("Sulfonamide", product)

                if is_sulfonamide_reaction or has_sulfonamide:
                    # Verify sulfonamide is formed
                    if has_sulfonamide:
                        sulfonamide_formed = True
                        sulfonamide_depth = depth
                        print(f"Sulfonamide formation detected at depth {depth}")

                        # Check if product has both sulfonamide and heterocycle
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check for heterocycles in the product
                            for heterocycle in heterocycles:
                                if checker.check_ring(heterocycle, product):
                                    # If both sulfonamide and heterocycle are present in the product,
                                    # they are likely connected (especially in a synthetic route)
                                    heterocycle_connected_to_sulfonamide = True
                                    print(
                                        f"Heterocycle {heterocycle} detected in product with sulfonamide"
                                    )
                                    break

                # If we haven't confirmed connection yet, check reactants for heterocycles
                if sulfonamide_formed and not heterocycle_connected_to_sulfonamide:
                    reactants = reactants_str.split(".")
                    for reactant in reactants:
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant):
                                # If a reactant has a heterocycle and we're forming a sulfonamide,
                                # the heterocycle is likely involved in the reaction
                                heterocycle_connected_to_sulfonamide = True
                                print(f"Heterocycle {heterocycle} detected in reactant: {reactant}")
                                break
                        if heterocycle_connected_to_sulfonamide:
                            break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if sulfonamide is formed late in the synthesis (depth ≤ 2) and heterocycles are present
    result = sulfonamide_formed and sulfonamide_depth <= 2 and heterocycle_connected_to_sulfonamide
    print(f"Late-stage sulfonamide formation strategy detected: {result}")
    return result
