#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis involves multiple C-N bond formations.
    """
    cn_bond_count = 0

    def dfs_traverse(node):
        nonlocal cn_bond_count

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]

            try:
                # Split reaction SMILES to analyze reactants and product
                parts = rsmi.split(">")
                reactants = parts[0].split(".")
                product = parts[-1]  # Use parts[-1] to safely get the product

                # Check if the reaction involves C-N bond formation
                is_cn_formation = False

                # Check named C-N bond formation reactions
                named_reactions = [
                    # SNAr and related reactions
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Ullmann-Goldberg Substitution amine",
                    "Goldberg coupling",
                    "Goldberg coupling aryl amine-aryl chloride",
                    "Goldberg coupling aryl amide-aryl chloride",
                    "Buchwald-Hartwig",
                    "N-arylation_heterocycles",
                    # Amide formation reactions
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
                    # Reductive amination
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "reductive amination",
                    # Urea formation
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                    "Urea synthesis via isocyanate and diazo",
                    "Urea synthesis via isocyanate and sulfonamide",
                    "urea",
                    # Other C-N bond formations
                    "Alkylation of amines",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "aza-Michael addition aromatic",
                    "aza-Michael addition secondary",
                    "aza-Michael addition primary",
                    "Aminolysis of esters",
                    "Ring opening of epoxide with amine",
                    "Intramolecular amination (heterocycle formation)",
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                ]

                for reaction_name in named_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_cn_formation = True
                        print(f"Found C-N bond formation via {reaction_name}: {rsmi}")
                        break

                # If not a named reaction, check for general C-N bond formation patterns
                if not is_cn_formation:
                    # Check for amine/amide/nitrogen-containing reactants
                    nitrogen_reactants = [
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Aniline",
                        "Azide",
                        "Hydrazine",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                    ]

                    has_nitrogen_reactant = False
                    for reactant in reactants:
                        for fg in nitrogen_reactants:
                            if checker.check_fg(fg, reactant):
                                has_nitrogen_reactant = True
                                break
                        if has_nitrogen_reactant:
                            break

                    # Check for carbon-containing electrophiles
                    carbon_electrophiles = [
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Acyl halide",
                        "Carboxylic acid",
                        "Ester",
                        "Aldehyde",
                        "Ketone",
                        "Nitrile",
                    ]

                    has_carbon_electrophile = False
                    for reactant in reactants:
                        for fg in carbon_electrophiles:
                            if checker.check_fg(fg, reactant):
                                has_carbon_electrophile = True
                                break
                        if has_carbon_electrophile:
                            break

                    # Check if product has new nitrogen-carbon bonds
                    # This is a simplified check - in a real implementation, we would need to analyze
                    # the atom mapping to confirm new C-N bonds are formed
                    nitrogen_carbon_products = [
                        "Aniline",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                        "Urea",
                        "Thiourea",
                        "Carbamate",
                    ]

                    has_nitrogen_carbon_product = False
                    for fg in nitrogen_carbon_products:
                        if checker.check_fg(fg, product):
                            has_nitrogen_carbon_product = True
                            break

                    if (
                        has_nitrogen_reactant
                        and has_carbon_electrophile
                        and has_nitrogen_carbon_product
                    ):
                        # This is a potential C-N bond formation
                        # In a more sophisticated implementation, we would check atom mappings
                        is_cn_formation = True
                        print(f"Found potential C-N bond formation: {rsmi}")

                if is_cn_formation:
                    cn_bond_count += 1
            except Exception as e:
                print(f"Error processing reaction SMILES: {rsmi}, Error: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total C-N bond formations found: {cn_bond_count}")
    return cn_bond_count >= 2  # Return True if at least 2 C-N bonds are formed
