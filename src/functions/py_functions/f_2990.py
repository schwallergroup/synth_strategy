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
    This function detects a synthetic strategy involving multiple C-N bond formations.
    """
    # Track C-N bond formations
    cn_bond_formations = 0

    # List of reaction types that form C-N bonds
    cn_bond_forming_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Aminolysis of esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with primary amine to imide",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Schotten-Baumann to ester",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acylation of secondary amines with anhydrides",
        "aza-Michael addition aromatic",
        "aza-Michael addition secondary",
        "aza-Michael addition primary",
        "Displacement of ethoxy group by primary amine",
        "Displacement of ethoxy group by secondary amine",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Sulfamoylarylamides from carboxylic acids and amines",
        "Carboxyl benzyl deprotection",
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "{reductive amination}",
    ]

    # Nitrogen-containing functional groups to check
    n_containing_fgs = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Aniline",
        "Azide",
        "Nitrile",
        "Nitro group",
        "Sulfonamide",
        "Urea",
        "Thiourea",
        "Isocyanate",
        "Substituted imine",
        "Unsubstituted imine",
    ]

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Examining reaction: {rsmi}")

                    # Extract reactants and product
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    # Check if this reaction forms a C-N bond using the checker
                    reaction_detected = False
                    for reaction_type in cn_bond_forming_reactions:
                        try:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(
                                    f"Detected C-N bond formation via {reaction_type}"
                                )
                                cn_bond_formations += 1
                                reaction_detected = True
                                break
                        except Exception as e:
                            print(f"Error checking reaction type {reaction_type}: {e}")

                    # If no reaction type was detected, check for N-containing functional group changes
                    if not reaction_detected:
                        # Check if product has N-containing functional groups not present in reactants
                        product_fgs = set()
                        reactant_fgs = set()

                        # Check functional groups in product
                        for fg in n_containing_fgs:
                            try:
                                if checker.check_fg(fg, product_part):
                                    product_fgs.add(fg)
                            except Exception as e:
                                print(f"Error checking FG {fg} in product: {e}")

                        # Check functional groups in reactants
                        for reactant in reactants:
                            for fg in n_containing_fgs:
                                try:
                                    if checker.check_fg(fg, reactant):
                                        reactant_fgs.add(fg)
                                except Exception as e:
                                    print(f"Error checking FG {fg} in reactant: {e}")

                        # If product has new N-containing functional groups, it might be a C-N bond formation
                        new_fgs = product_fgs - reactant_fgs
                        if new_fgs:
                            print(
                                f"Detected potential C-N bond formation: new functional groups {new_fgs}"
                            )
                            cn_bond_formations += 1
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    try:
        dfs_traverse(route)
        print(f"Total C-N bond formations detected: {cn_bond_formations}")
    except Exception as e:
        print(f"Error traversing route: {e}")

    # Check if we have multiple C-N bond formations
    if cn_bond_formations >= 2:
        print(
            f"Strategy detected: Multiple C-N bond formations ({cn_bond_formations} instances)"
        )
        return True

    return False
