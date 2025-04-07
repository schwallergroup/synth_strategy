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
    This function detects a synthesis route with multiple C-N bond formations.
    """
    cn_bond_formation_count = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formation_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for C-N bond formation reactions
                cn_bond_formation_reactions = [
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
                    "Schotten-Baumann_amide",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "N-arylation_heterocycles",
                    "Buchwald-Hartwig",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "reductive amination",
                    "Alkylation of amines",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                    "urea",
                    "Goldberg coupling aryl amine-aryl chloride",
                    "Goldberg coupling aryl amide-aryl chloride",
                    "Goldberg coupling",
                    "Ullmann-Goldberg Substitution amine",
                    "Aminolysis of esters",
                    "aza-Michael addition aromatic",
                    "aza-Michael addition secondary",
                    "aza-Michael addition primary",
                    "Ring opening of epoxide with amine",
                ]

                for reaction_type in cn_bond_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        cn_bond_formation_count += 1
                        print(f"Found C-N bond formation ({reaction_type}): {rsmi}")
                        break  # Count each reaction only once

                # Additional check for reactions not in the predefined list
                if not any(
                    checker.check_reaction(reaction_type, rsmi)
                    for reaction_type in cn_bond_formation_reactions
                ):
                    # Check for reactants with nitrogen-containing groups
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Check if reactants contain nitrogen groups
                    reactant_mols = reactants_part.split(".")
                    product_mol = product_part

                    n_groups_in_reactants = False
                    for reactant in reactant_mols:
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Azide", reactant)
                            or checker.check_fg("Hydrazine", reactant)
                        ):
                            n_groups_in_reactants = True
                            break

                    # Check if product has new C-N bonds not present in reactants
                    if n_groups_in_reactants:
                        # Check for amide formation
                        if (
                            checker.check_fg("Primary amide", product_mol)
                            or checker.check_fg("Secondary amide", product_mol)
                            or checker.check_fg("Tertiary amide", product_mol)
                        ):
                            cn_bond_formation_count += 1
                            print(f"Found C-N bond formation (amide formation): {rsmi}")
                        # Check for other C-N bond formations
                        elif checker.check_fg(
                            "Substituted imine", product_mol
                        ) or checker.check_fg("Unsubstituted imine", product_mol):
                            cn_bond_formation_count += 1
                            print(f"Found C-N bond formation (imine formation): {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 C-N bond formations are detected
    strategy_present = cn_bond_formation_count >= 2
    print(f"Found {cn_bond_formation_count} C-N bond formations")
    print(f"Strategy detection result: {strategy_present}")
    return strategy_present
