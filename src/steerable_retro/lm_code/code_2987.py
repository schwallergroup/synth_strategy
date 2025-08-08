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
    This function detects if the final step in the synthesis is an amide formation.
    """
    final_step_is_amide_formation = False

    def dfs_traverse_from_mol(mol_node):
        """Helper function to find the final reaction step when starting from a molecule node"""
        nonlocal final_step_is_amide_formation

        # Check if this molecule has any children (reactions)
        if "children" in mol_node and mol_node["children"]:
            for child in mol_node["children"]:
                if child["type"] == "reaction":
                    # This is the final reaction step
                    analyze_reaction(child, is_final_step=True)
                    # Only traverse further if we haven't found an amide formation
                    if not final_step_is_amide_formation:
                        for grandchild in child.get("children", []):
                            dfs_traverse(grandchild, 1)

    def analyze_reaction(reaction_node, is_final_step=False):
        """Analyze a reaction node to determine if it's an amide formation"""
        nonlocal final_step_is_amide_formation

        if "rsmi" in reaction_node["metadata"]:
            rsmi = reaction_node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing {'final' if is_final_step else 'intermediate'} reaction: {rsmi}")

            # Check if this is an amide formation reaction directly
            amide_formation_reactions = [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Acyl chloride with ammonia to amide",
                "Ester with ammonia to amide",
                "Schotten-Baumann to ester",
                "Schotten-Baumann_amide",
                "Acylation of secondary amines with anhydrides",
                "Acylation of primary amines",
                "Acylation of secondary amines",
            ]

            for reaction_type in amide_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(
                        f"Detected {'late-stage' if is_final_step else 'intermediate'} amide formation: {reaction_type}"
                    )
                    if is_final_step:
                        final_step_is_amide_formation = True
                    return

            print("No direct amide formation reaction type matched, checking functional groups")

            # If no specific reaction type matched, check for functional group changes
            # Check if product contains amide group
            has_amide_product = False
            for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                if checker.check_fg(fg, product):
                    has_amide_product = True
                    print(f"Found {fg} in product")

            if not has_amide_product:
                print("No amide in product, not an amide formation")
                return  # No amide in product, not an amide formation

            # Check if reactants contain carboxylic acid/derivatives and amine
            has_carboxylic_acid = False
            has_acyl_chloride = False
            has_ester = False
            has_anhydride = False
            has_amine = False

            for reactant in reactants:
                if checker.check_fg("Carboxylic acid", reactant):
                    has_carboxylic_acid = True
                    print("Found carboxylic acid in reactants")
                if checker.check_fg("Acyl halide", reactant):
                    has_acyl_chloride = True
                    print("Found acyl halide in reactants")
                if checker.check_fg("Ester", reactant):
                    has_ester = True
                    print("Found ester in reactants")
                if checker.check_fg("Anhydride", reactant):
                    has_anhydride = True
                    print("Found anhydride in reactants")
                if checker.check_fg("Primary amine", reactant):
                    has_amine = True
                    print("Found primary amine in reactants")
                if checker.check_fg("Secondary amine", reactant):
                    has_amine = True
                    print("Found secondary amine in reactants")
                if checker.check_fg("Tertiary amine", reactant):
                    has_amine = True
                    print("Found tertiary amine in reactants")

            # If product has amide and reactants have acid/derivatives and amine, it's amide formation
            if (
                has_amide_product
                and has_amine
                and (has_carboxylic_acid or has_acyl_chloride or has_ester or has_anhydride)
            ):
                print(
                    f"Detected {'late-stage' if is_final_step else 'intermediate'} amide formation from functional group analysis"
                )
                if is_final_step:
                    final_step_is_amide_formation = True

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_formation

        # If root is a molecule, we need special handling
        if node["type"] == "mol" and depth == 0:
            print("Root node is a molecule, finding final reaction step")
            dfs_traverse_from_mol(node)
            return

        if node["type"] == "reaction" and depth == 0:  # Final step (depth 0)
            analyze_reaction(node, is_final_step=True)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return final_step_is_amide_formation
