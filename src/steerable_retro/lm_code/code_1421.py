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
    Detects aldehyde reduction to alcohol as part of the synthetic strategy.
    """
    found_aldehyde_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_aldehyde_reduction

        if found_aldehyde_reduction:
            return  # Early return if already found

        if node["type"] == "reaction":
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Examining reaction at depth {depth}: {rsmi}")

                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")
                    products = product_part.split(".")

                    # Check for reduction reaction (forward direction)
                    if checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    ):
                        print(
                            f"Matched 'Reduction of aldehydes and ketones to alcohols' reaction type"
                        )

                        # Check for aldehyde in reactants and primary alcohol in products
                        for reactant in reactants:
                            if checker.check_fg("Aldehyde", reactant):
                                print(f"Found aldehyde in reactant: {reactant}")
                                for product in products:
                                    if checker.check_fg("Primary alcohol", product):
                                        print(f"Found primary alcohol in product: {product}")
                                        found_aldehyde_reduction = True
                                        return

                    # Check for oxidation reaction (reverse direction)
                    elif checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    ):
                        print(
                            f"Matched 'Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones' reaction type"
                        )

                        # In retrosynthesis, this would be a reduction when traversing backward
                        for reactant in reactants:
                            if checker.check_fg("Primary alcohol", reactant):
                                print(f"Found primary alcohol in reactant: {reactant}")
                                for product in products:
                                    if checker.check_fg("Aldehyde", product):
                                        print(f"Found aldehyde in product: {product}")
                                        found_aldehyde_reduction = True
                                        return

                    # Check for uncategorized reactions that might involve aldehyde reduction
                    else:
                        # Check for aldehyde in reactants and primary alcohol in products
                        for reactant in reactants:
                            if checker.check_fg("Aldehyde", reactant):
                                print(f"Found aldehyde in reactant: {reactant}")
                                for product in products:
                                    if checker.check_fg(
                                        "Primary alcohol", product
                                    ) and not checker.check_fg("Aldehyde", product):
                                        print(
                                            f"Found potential aldehyde reduction to alcohol (uncategorized): {rsmi}"
                                        )

                                        # Exclude other types of reductions
                                        if not (
                                            checker.check_reaction(
                                                "Reduction of carboxylic acid to primary alcohol",
                                                rsmi,
                                            )
                                            or checker.check_reaction(
                                                "Reduction of ester to primary alcohol", rsmi
                                            )
                                            or checker.check_reaction(
                                                "Reduction of nitrile to amine", rsmi
                                            )
                                        ):
                                            print(
                                                f"Confirmed aldehyde reduction to alcohol: {rsmi}"
                                            )
                                            found_aldehyde_reduction = True
                                            return

                        # Check for primary alcohol in reactants and aldehyde in products (reverse direction)
                        for reactant in reactants:
                            if checker.check_fg("Primary alcohol", reactant):
                                print(f"Found primary alcohol in reactant: {reactant}")
                                for product in products:
                                    if checker.check_fg(
                                        "Aldehyde", product
                                    ) and not checker.check_fg("Primary alcohol", product):
                                        print(
                                            f"Found potential alcohol oxidation to aldehyde (uncategorized): {rsmi}"
                                        )
                                        found_aldehyde_reduction = True
                                        return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children nodes (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_aldehyde_reduction
