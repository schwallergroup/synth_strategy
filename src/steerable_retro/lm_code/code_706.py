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
    This function detects a linear synthesis strategy with multiple C-N bond formations.
    """
    cn_bond_formations = 0
    reaction_count = 0

    # List of reaction types that form C-N bonds
    cn_bond_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "N-methylation",
        "Aminolysis of esters",
        "Schotten-Baumann_amide",
        "Buchwald-Hartwig",
        "reductive amination",
    ]

    # C-N containing functional groups to check
    cn_fgs = [
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Urea",
        "Sulfonamide",
        "Aniline",
        "Amidinium",
        "Carbamic ester",
        "Carbamic acid",
        "Cyanamide",
    ]

    def is_linear_path(node):
        """Check if a node is part of a linear path (has at most one child)"""
        if node["type"] == "mol" and not node.get("in_stock", False):
            # Non-terminal molecule nodes should have exactly one child (reaction)
            return len(node.get("children", [])) == 1
        elif node["type"] == "reaction":
            # Reaction nodes can have multiple children (reactants)
            # But at least one of those children should continue the linear path
            linear_children = 0
            for child in node.get("children", []):
                if is_linear_path(child):
                    linear_children += 1
            return linear_children >= 1
        return True  # Terminal nodes (in_stock molecules) are considered part of a linear path

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and product
            try:
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]
                reactants = reactants_str.split(".")

                # Method 1: Check if this reaction is a known C-N bond forming reaction
                for rxn_type in cn_bond_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(
                            f"Found C-N bond formation via reaction type: {rxn_type} at depth {depth}"
                        )
                        cn_bond_formations += 1
                        break

                # Method 2: Check for appearance of C-N functional groups
                else:
                    # Check if product has C-N functional groups not present in reactants
                    for fg in cn_fgs:
                        if checker.check_fg(fg, product_str):
                            # Check if this FG was not in any reactant
                            if not any(checker.check_fg(fg, r) for r in reactants):
                                print(
                                    f"Found C-N bond formation via new functional group: {fg} at depth {depth}"
                                )
                                cn_bond_formations += 1
                                break
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Check if the route is linear
    if not is_linear_path(route):
        print("Route is not linear")
        return False

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Total reactions: {reaction_count}, C-N bond formations: {cn_bond_formations}")

    # Return True if we have multiple C-N bond formations in a linear synthesis
    # (number of C-N formations is at least 2 and close to the number of reactions)
    return cn_bond_formations >= 2 and cn_bond_formations >= reaction_count * 0.5
