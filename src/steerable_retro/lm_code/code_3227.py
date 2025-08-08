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
    This function detects nitro group reduction to amine during synthesis.
    """
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Direct check for nitro reduction reaction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Found nitro reduction reaction at depth {depth}")
                nitro_reduction_found = True
                return

            # Check for nitro groups in reactants
            for reactant in reactants:
                if checker.check_fg("Nitro group", reactant):
                    print(f"Found nitro group in reactant: {reactant}")

                    # Check if the nitro group is absent in the product
                    if not checker.check_fg("Nitro group", product):
                        print(f"Nitro group is absent in product: {product}")

                        # Check for various amine types in the product
                        if (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                        ):

                            print(f"Found amine in product: {product}")

                            # Additional check for reaction types that might involve nitro reduction
                            if (
                                checker.check_reaction("Reduction of nitro groups to amines", rsmi)
                                or "reduction" in rsmi.lower()
                                or "nitro" in rsmi.lower()
                                and "amine" in rsmi.lower()
                            ):
                                print(f"Found nitro group reduction to amine at depth {depth}")
                                nitro_reduction_found = True
                                return

                    # Even if nitro group is still present, check if any nitro group was reduced
                    # This handles partial reduction in molecules with multiple nitro groups
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Count nitro groups in reactant and product
                        reactant_nitro_count = len(
                            checker.get_fg_atom_indices("Nitro group", reactant)
                        )
                        product_nitro_count = len(
                            checker.get_fg_atom_indices("Nitro group", product)
                        )

                        # If there are fewer nitro groups in the product, and new amines appeared
                        if product_nitro_count < reactant_nitro_count and (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                        ):

                            print(f"Found partial nitro group reduction at depth {depth}")
                            print(
                                f"Reactant nitro count: {reactant_nitro_count}, Product nitro count: {product_nitro_count}"
                            )
                            nitro_reduction_found = True
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitro_reduction_found
