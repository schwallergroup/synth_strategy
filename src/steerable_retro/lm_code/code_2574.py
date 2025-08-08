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
    Detects a synthetic strategy involving late-stage allylation,
    particularly focusing on O-alkylation with allyl groups.
    """
    # Initialize tracking variables
    has_late_allylation = False

    # Add depth information to the route
    def add_depth(node, current_depth=0):
        node["depth"] = current_depth
        for child in node.get("children", []):
            add_depth(child, current_depth + 1)

    add_depth(route)

    def dfs_traverse(node):
        nonlocal has_late_allylation

        if node["type"] == "reaction":
            # Consider reactions in the late stage (depth <= 2)
            if node.get("depth", 0) <= 2:
                print(f"Examining reaction at depth {node.get('depth', 0)}")

                # Extract reactants and products
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants_smiles}")
                    print(f"Product: {product_smiles}")

                    # Check if this is a reaction that could involve allylation
                    allylation_reaction = (
                        checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction("S-alkylation of thiols", rsmi)
                        or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                        or checker.check_reaction(
                            "Williamson Ether Synthesis (intra to epoxy)", rsmi
                        )
                    )

                    # First check: Is the product an allyl ether/thioether?
                    if checker.check_fg("Allyl", product_smiles) and (
                        checker.check_fg("Ether", product_smiles)
                        or checker.check_fg("Monosulfide", product_smiles)
                    ):

                        print("Product contains allyl ether/thioether")

                        # Second check: Is this a known allylation reaction or does it follow the pattern?
                        if allylation_reaction:
                            print("Detected known allylation reaction")
                            has_late_allylation = True
                        else:
                            # Check if any reactant has OH/SH group (potential allylation substrate)
                            has_oh_sh_reactant = False
                            has_allyl_reactant = False

                            for reactant in reactants_smiles:
                                # Check for OH/SH containing reactant
                                if (
                                    checker.check_fg("Phenol", reactant)
                                    or checker.check_fg("Primary alcohol", reactant)
                                    or checker.check_fg("Secondary alcohol", reactant)
                                    or checker.check_fg("Tertiary alcohol", reactant)
                                    or checker.check_fg("Aliphatic thiol", reactant)
                                    or checker.check_fg("Aromatic thiol", reactant)
                                ):
                                    has_oh_sh_reactant = True
                                    print(f"Found reactant with OH/SH group: {reactant}")

                                # Check for allyl-containing reactant
                                if checker.check_fg("Allyl", reactant):
                                    has_allyl_reactant = True
                                    print(f"Found allyl group in reactant: {reactant}")

                            # If we have both an OH/SH reactant and an allyl reactant, this is likely allylation
                            if has_oh_sh_reactant and has_allyl_reactant:
                                print(
                                    f"Detected pattern-based late-stage allylation at depth {node.get('depth', 0)}"
                                )
                                has_late_allylation = True

                            # Special case: If the allyl group is already present in an OH/SH reactant,
                            # and preserved in the product, this could be a rearrangement or protection
                            elif not has_late_allylation:
                                for reactant in reactants_smiles:
                                    if checker.check_fg("Allyl", reactant) and (
                                        checker.check_fg("Ether", reactant)
                                        or checker.check_fg("Monosulfide", reactant)
                                    ):
                                        print(
                                            f"Found allyl ether/thioether in reactant: {reactant}"
                                        )
                                        # This is a reaction involving an allyl ether/thioether, which counts as late-stage allylation
                                        has_late_allylation = True
                                        print(
                                            f"Detected late-stage allyl ether/thioether transformation at depth {node.get('depth', 0)}"
                                        )
                                        break

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage allylation strategy detected: {has_late_allylation}")
    return has_late_allylation
