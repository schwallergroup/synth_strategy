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
    This function detects a linear synthesis strategy with multiple C-N bond formations,
    including amine alkylation, guanidine/urea formation, and amide coupling.
    """
    cn_bond_formations = 0
    has_amide_coupling = False
    has_guanidine_formation = False
    has_amine_alkylation = False
    is_linear = True

    # List of reaction types that form C-N bonds
    cn_bond_forming_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Alkylation of amines",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Buchwald-Hartwig",
        "Aminolysis of esters",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Addition of primary amines to aldehydes/thiocarbonyls",
        "Addition of primary amines to ketones/thiocarbonyls",
        "Addition of secondary amines to ketones/thiocarbonyls",
        "Addition of secondary amines to aldehydes/thiocarbonyls",
        "Mignonac reaction",
        "Reductive amination",
        "aza-Michael addition aromatic",
        "aza-Michael addition secondary",
        "aza-Michael addition primary",
        "Ugi reaction",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Paal-Knorr pyrrole synthesis",
        "Pictet-Spengler",
        "Urea",
        "Thiourea",
    ]

    # Reaction types by category
    amide_coupling_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
        "Acylation of primary amines",
        "Acylation of secondary amines",
    ]

    guanidine_urea_reactions = [
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
        "Urea",
        "Thiourea",
    ]

    amine_alkylation_reactions = [
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Alkylation of amines",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Reductive amination",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Buchwald-Hartwig",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
    ]

    # Additional functional groups to check
    amine_fgs = ["Primary amine", "Secondary amine", "Tertiary amine", "Aniline"]
    amide_fgs = ["Primary amide", "Secondary amide", "Tertiary amide"]
    urea_fgs = ["Urea", "Thiourea"]

    def has_cn_bond_formation(rsmi, reactants, product):
        """Check if a reaction forms a C-N bond by comparing reactants and product"""
        # First check if it's a known C-N bond forming reaction
        for reaction_type in cn_bond_forming_reactions:
            if checker.check_reaction(reaction_type, rsmi):
                print(f"Detected C-N bond formation reaction: {reaction_type}")
                return True, reaction_type

        # If not a known reaction, try to detect C-N bond formation by checking
        # for amine/amide/urea functional groups in the product that weren't in reactants
        try:
            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return False, None

            # Check for new amine/amide/urea groups in product
            for fg in amine_fgs + amide_fgs + urea_fgs:
                if checker.check_fg(fg, product):
                    # Check if this FG was already in reactants
                    fg_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg(fg, reactant):
                            fg_in_reactants = True
                            break

                    if not fg_in_reactants:
                        print(f"Detected new {fg} functional group in product")
                        # Determine reaction category
                        if fg in amide_fgs:
                            return True, "Amide formation"
                        elif fg in urea_fgs:
                            return True, "Urea/Guanidine formation"
                        elif fg in amine_fgs:
                            return True, "Amine alkylation"

            return False, None
        except Exception as e:
            print(f"Error checking for C-N bond formation: {e}")
            return False, None

    def is_synthesis_linear(reactants, product):
        """Determine if a synthesis step is linear by comparing reactants and product"""
        try:
            # If there's only one reactant, it's linear
            if len(reactants) <= 1:
                return True

            # If there are two reactants, check if one is significantly smaller (reagent)
            if len(reactants) == 2:
                mol1 = Chem.MolFromSmiles(reactants[0])
                mol2 = Chem.MolFromSmiles(reactants[1])

                if not mol1 or not mol2:
                    return True  # Default to True if parsing fails

                # Count heavy atoms
                count1 = mol1.GetNumHeavyAtoms()
                count2 = mol2.GetNumHeavyAtoms()

                # If one reactant is much smaller, consider it a reagent
                if count1 > 2 * count2 or count2 > 2 * count1:
                    return True

                # If both reactants are reasonably sized, still consider it linear
                # as most organic reactions with two main components are still linear
                if count1 >= 5 and count2 >= 5:
                    return True

            # For more than 2 reactants, most organic reactions are still considered linear
            # in the context of synthetic planning
            return True

        except Exception as e:
            print(f"Error checking if synthesis is linear: {e}")
            return True  # Default to True if there's an error

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations, has_amide_coupling, has_guanidine_formation, has_amine_alkylation, is_linear

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                print(f"\nAnalyzing reaction at depth {depth}: {rsmi}")

                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if synthesis is linear
                if not is_synthesis_linear(reactants, product_part):
                    is_linear = False
                    print(f"Not a linear synthesis at depth {depth}: {rsmi}")

                # Check for C-N bond formation
                formed_cn_bond, reaction_category = has_cn_bond_formation(
                    rsmi, reactants, product_part
                )

                if formed_cn_bond:
                    cn_bond_formations += 1
                    print(f"C-N bond formation detected at depth {depth}")

                    # Categorize the reaction
                    if (
                        reaction_category in amide_coupling_reactions
                        or reaction_category == "Amide formation"
                    ):
                        has_amide_coupling = True
                        print(f"Categorized as amide coupling at depth {depth}")
                    elif (
                        reaction_category in guanidine_urea_reactions
                        or reaction_category == "Urea/Guanidine formation"
                    ):
                        has_guanidine_formation = True
                        print(f"Categorized as guanidine/urea formation at depth {depth}")
                    elif (
                        reaction_category in amine_alkylation_reactions
                        or reaction_category == "Amine alkylation"
                    ):
                        has_amine_alkylation = True
                        print(f"Categorized as amine alkylation at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction node at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        cn_bond_formations >= 2
        and is_linear
        and (has_amide_coupling or has_guanidine_formation or has_amine_alkylation)
    )

    print(f"\nSummary:")
    print(f"Linear synthesis with multiple C-N bonds: {strategy_present}")
    print(f"C-N bond formations: {cn_bond_formations}")
    print(f"Is linear: {is_linear}")
    print(f"Has amide coupling: {has_amide_coupling}")
    print(f"Has guanidine formation: {has_guanidine_formation}")
    print(f"Has amine alkylation: {has_amine_alkylation}")

    return strategy_present
