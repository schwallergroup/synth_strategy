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
    This function detects the strategy of introducing a 4-nitrobenzenesulfonyl group.
    """
    has_nitro_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_sulfonamide

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a sulfonamide formation reaction
                is_sulfonamide_reaction = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                # Also check for general sulfonamide formation
                if not is_sulfonamide_reaction:
                    is_sulfonamide_reaction = checker.check_fg("Sulfonamide", product) and any(
                        checker.check_fg("Sulfonyl halide", reactant) for reactant in reactants
                    )

                if is_sulfonamide_reaction:
                    # Check if product contains sulfonamide group
                    has_sulfonamide = checker.check_fg("Sulfonamide", product)

                    # Check if one of the reactants contains a nitro group on an aromatic ring
                    nitro_sulfonyl_reactant = None
                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                            "Sulfonyl halide", reactant
                        ):
                            nitro_sulfonyl_reactant = reactant
                            break

                    has_nitro_reactant = nitro_sulfonyl_reactant is not None

                    # Check if the product has both nitro and sulfonamide groups
                    has_nitro_product = checker.check_fg("Nitro group", product)

                    # If we have a sulfonamide formation with a nitro group preserved from reactant to product
                    if has_sulfonamide and has_nitro_reactant and has_nitro_product:
                        # Check if the product has a benzene ring
                        has_benzene = checker.check_ring("benzene", product)

                        if has_benzene:
                            print(f"Found nitrobenzenesulfonamide formation at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            print(f"Nitro-sulfonyl reactant: {nitro_sulfonyl_reactant}")
                            print(f"Product: {product}")
                            has_nitro_sulfonamide = True

                # Also check for direct use of 4-nitrobenzenesulfonamide in any reaction
                for reactant in reactants:
                    if (
                        checker.check_fg("Nitro group", reactant)
                        and checker.check_fg("Sulfonamide", reactant)
                        and checker.check_ring("benzene", reactant)
                    ):
                        print(f"Found direct use of nitrobenzenesulfonamide at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Nitro-sulfonamide reactant: {reactant}")
                        has_nitro_sulfonamide = True
                        break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Also check if the molecule itself is a nitrobenzenesulfonamide
        elif node["type"] == "mol" and not node.get("in_stock", False):
            mol_smiles = node["smiles"]
            if (
                checker.check_fg("Nitro group", mol_smiles)
                and checker.check_fg("Sulfonamide", mol_smiles)
                and checker.check_ring("benzene", mol_smiles)
            ):
                print(f"Found nitrobenzenesulfonamide molecule at depth {depth}")
                print(f"Molecule SMILES: {mol_smiles}")
                has_nitro_sulfonamide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    print(f"Nitro sulfonamide strategy detected: {has_nitro_sulfonamide}")
    return has_nitro_sulfonamide
