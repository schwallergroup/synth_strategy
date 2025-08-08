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
    Detects if the synthetic route involves attaching an indole moiety
    in the final step via C-O bond formation.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        print(f"Examining node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is the final reaction (depth 0 or 1)
            node_depth = node["metadata"].get("depth", depth)
            print(f"Reaction at depth {node_depth}: {rsmi}")

            if node_depth <= 1:  # Consider both depth 0 and 1 as potentially final reactions
                print(f"Examining potential final reaction: {rsmi}")

                # Check for indole in reactants
                indole_in_reactants = False
                indole_reactant = None

                for reactant in reactants_smiles:
                    try:
                        if checker.check_ring("indole", reactant):
                            print(f"Found indole in reactant: {reactant}")
                            indole_in_reactants = True
                            indole_reactant = reactant
                            break
                    except Exception as e:
                        print(f"Error checking for indole in reactant: {e}")

                # Check if product contains indole
                indole_in_product = False
                try:
                    indole_in_product = checker.check_ring("indole", product_smiles)
                    print(f"Indole in product: {indole_in_product}")
                except Exception as e:
                    print(f"Error checking for indole in product: {e}")

                if indole_in_reactants and indole_reactant and indole_in_product:
                    # Check if this is a C-O bond forming reaction
                    co_bond_reaction = False

                    # Check for common C-O bond forming reactions
                    if (
                        checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                        or checker.check_reaction("Chan-Lam etherification", rsmi)
                        or checker.check_reaction(
                            "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                        )
                        or checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                        or checker.check_reaction("Alcohol to ether", rsmi)
                        or checker.check_reaction("{Williamson ether}", rsmi)
                        or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                        or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                        or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    ):
                        print(f"Found C-O bond forming reaction: {rsmi}")
                        co_bond_reaction = True

                    # If no specific reaction type is detected, check if this is any reaction forming a C-O bond
                    if not co_bond_reaction:
                        # Check if the reaction involves a nucleophilic substitution with oxygen
                        for reactant in reactants_smiles:
                            if (
                                checker.check_fg("Phenol", reactant)
                                or checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Tertiary alcohol", reactant)
                                or checker.check_fg("Aromatic alcohol", reactant)
                            ):
                                for other_reactant in reactants_smiles:
                                    if other_reactant != reactant and (
                                        checker.check_fg("Aromatic halide", other_reactant)
                                        or checker.check_fg("Primary halide", other_reactant)
                                    ):
                                        if checker.check_fg("Ether", product_smiles):
                                            print("Detected general C-O bond formation")
                                            co_bond_reaction = True
                                            break

                    # Check for ether formation in the product
                    if co_bond_reaction or checker.check_fg("Ether", product_smiles):
                        # Check if the indole reactant contains a nucleophilic oxygen group
                        indole_has_oxygen_nucleophile = (
                            checker.check_fg("Phenol", indole_reactant)
                            or checker.check_fg("Primary alcohol", indole_reactant)
                            or checker.check_fg("Secondary alcohol", indole_reactant)
                            or checker.check_fg("Tertiary alcohol", indole_reactant)
                            or checker.check_fg("Aromatic alcohol", indole_reactant)
                        )

                        if indole_has_oxygen_nucleophile:
                            print("Indole reactant contains oxygen nucleophile")
                            # Check if product has ether group
                            if checker.check_fg("Ether", product_smiles):
                                print("Product contains ether group - pattern found!")
                                found_pattern = True
                        else:
                            # Special case for phenols which can act as nucleophiles in SNAr reactions
                            if checker.check_fg("Phenol", indole_reactant):
                                print("Indole reactant contains phenol group")
                                if checker.check_fg("Ether", product_smiles):
                                    print("Product contains ether group - pattern found!")
                                    found_pattern = True

                            # Check if indole has a leaving group (electrophile)
                            indole_has_leaving_group = (
                                checker.check_fg("Primary halide", indole_reactant)
                                or checker.check_fg("Secondary halide", indole_reactant)
                                or checker.check_fg("Tertiary halide", indole_reactant)
                                or checker.check_fg("Aromatic halide", indole_reactant)
                                or checker.check_fg("Triflate", indole_reactant)
                                or checker.check_fg("Mesylate", indole_reactant)
                                or checker.check_fg("Tosylate", indole_reactant)
                            )

                            if indole_has_leaving_group:
                                print("Indole reactant contains leaving group")
                                # Check if another reactant has an oxygen nucleophile
                                for other_reactant in reactants_smiles:
                                    if other_reactant != indole_reactant and (
                                        checker.check_fg("Phenol", other_reactant)
                                        or checker.check_fg("Primary alcohol", other_reactant)
                                        or checker.check_fg("Secondary alcohol", other_reactant)
                                        or checker.check_fg("Tertiary alcohol", other_reactant)
                                        or checker.check_fg("Aromatic alcohol", other_reactant)
                                    ):
                                        print(
                                            f"Found oxygen nucleophile in other reactant: {other_reactant}"
                                        )
                                        # Check if product has ether group
                                        if checker.check_fg("Ether", product_smiles):
                                            print("Product contains ether group - pattern found!")
                                            found_pattern = True
                                            break
                            else:
                                # Check if any other reactant has both oxygen and leaving group
                                for other_reactant in reactants_smiles:
                                    if other_reactant != indole_reactant:
                                        # Check if this reactant has a leaving group
                                        other_has_leaving_group = (
                                            checker.check_fg("Primary halide", other_reactant)
                                            or checker.check_fg("Secondary halide", other_reactant)
                                            or checker.check_fg("Tertiary halide", other_reactant)
                                            or checker.check_fg("Aromatic halide", other_reactant)
                                            or checker.check_fg("Triflate", other_reactant)
                                            or checker.check_fg("Mesylate", other_reactant)
                                            or checker.check_fg("Tosylate", other_reactant)
                                        )

                                        if other_has_leaving_group:
                                            print(
                                                f"Found leaving group in other reactant: {other_reactant}"
                                            )
                                            # Check if product has ether group
                                            if checker.check_fg("Ether", product_smiles):
                                                print(
                                                    "Product contains ether group - pattern found!"
                                                )
                                                found_pattern = True
                                                break

        # Continue traversing
        for i, child in enumerate(node.get("children", [])):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Final result: {found_pattern}")
    return found_pattern
