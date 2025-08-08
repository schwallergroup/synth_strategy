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
    This function detects if the synthesis involves multiple C-O bond formations or cleavages.
    """
    c_o_bond_operations = 0

    def dfs_traverse(node):
        nonlocal c_o_bond_operations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Extract reactants and product
                try:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    # Check for C-O bond formation reactions
                    c_o_formation_reactions = [
                        "Williamson Ether Synthesis",
                        "Esterification of Carboxylic Acids",
                        "Alcohol protection with silyl ethers",
                        "O-alkylation of carboxylic acids with diazo compounds",
                        "Aldehyde or ketone acetalization",
                        "Diol acetalization",
                        "Oxidative esterification of primary alcohols",
                        "Alcohol to ether",
                        "Mitsunobu aryl ether",
                        "Mitsunobu esterification",
                    ]

                    # Check for C-O bond cleavage reactions
                    c_o_cleavage_reactions = [
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "Alcohol deprotection from silyl ethers",
                        "Acetal hydrolysis to diol",
                        "Acetal hydrolysis to aldehyde",
                        "Ketal hydrolysis to ketone",
                        "Cleavage of methoxy ethers to alcohols",
                        "Cleavage of alkoxy ethers to alcohols",
                        "Ether cleavage to primary alcohol",
                    ]

                    # Check if the reaction involves C-O bond formation or cleavage
                    for rxn_type in c_o_formation_reactions + c_o_cleavage_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"C-O bond manipulation detected: {rxn_type}")
                            c_o_bond_operations += 1
                            break

                    # If no specific reaction type was found, check for functional group changes
                    # that indicate C-O bond manipulation
                    if not any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in c_o_formation_reactions + c_o_cleavage_reactions
                    ):
                        # Check for ether formation
                        product_has_ether = checker.check_fg("Ether", product_part)
                        reactants_have_alcohol = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            for r in reactants
                        )

                        # Check for ester formation
                        product_has_ester = checker.check_fg("Ester", product_part)
                        reactants_have_carboxylic_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )
                        reactants_have_alcohol = reactants_have_alcohol or any(
                            checker.check_fg("Aromatic alcohol", r) for r in reactants
                        )

                        # Check for alcohol deprotection
                        product_has_alcohol = (
                            checker.check_fg("Primary alcohol", product_part)
                            or checker.check_fg("Secondary alcohol", product_part)
                            or checker.check_fg("Tertiary alcohol", product_part)
                            or checker.check_fg("Aromatic alcohol", product_part)
                        )

                        # Detect C-O bond formation
                        if (product_has_ether and reactants_have_alcohol) or (
                            product_has_ester
                            and reactants_have_carboxylic_acid
                            and reactants_have_alcohol
                        ):
                            print(f"C-O bond formation detected through functional group analysis")
                            c_o_bond_operations += 1

                        # Detect C-O bond cleavage
                        elif product_has_alcohol and (
                            any(checker.check_fg("Ether", r) for r in reactants)
                            or any(checker.check_fg("Ester", r) for r in reactants)
                        ):
                            print(f"C-O bond cleavage detected through functional group analysis")
                            c_o_bond_operations += 1

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total C-O bond operations detected: {c_o_bond_operations}")
    # Return True if at least 2 C-O bond operations are detected
    return c_o_bond_operations >= 2
