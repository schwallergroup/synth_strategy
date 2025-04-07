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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects if the route incorporates a protected amine (specifically Boc-protected)
    in a coupling reaction.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth <= 1)
            if depth <= 1:
                # Check if any reactant contains a Boc-protected amine
                boc_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Boc", reactant):
                        boc_in_reactants = True
                        print(f"Found Boc-protected amine in reactant: {reactant}")
                        break

                # Check if product contains Boc-protected amine
                boc_in_product = checker.check_fg("Boc", product)
                if boc_in_product:
                    print(f"Found Boc-protected amine in product: {product}")

                # If Boc is preserved through the reaction
                if boc_in_reactants and boc_in_product:
                    # Check if this is a coupling reaction
                    coupling_reactions = [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic esters",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Sonogashira acetylene_aryl halide",
                        "Sonogashira alkyne_aryl halide",
                        "Heck terminal vinyl",
                        "Negishi coupling",
                        "Stille reaction_aryl",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Schotten-Baumann to ester",
                        "N-alkylation of primary amines with alkyl halides",
                        "N-alkylation of secondary amines with alkyl halides",
                        "Alkylation of amines",
                    ]

                    is_coupling = False
                    for reaction_type in coupling_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_coupling = True
                            print(f"Found coupling reaction: {reaction_type}")
                            break

                    # Check for general N-alkylation if not found in specific types
                    if not is_coupling:
                        # Check if this is an N-alkylation reaction by looking for amine and halide
                        amine_in_reactants = False
                        halide_in_reactants = False

                        for reactant in reactants:
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                amine_in_reactants = True
                            if (
                                checker.check_fg("Primary halide", reactant)
                                or checker.check_fg("Secondary halide", reactant)
                                or checker.check_fg("Tertiary halide", reactant)
                                or checker.check_fg("Aromatic halide", reactant)
                            ):
                                halide_in_reactants = True

                        # If we have both amine and halide in reactants, it might be an N-alkylation
                        if amine_in_reactants and halide_in_reactants:
                            # Check if product has a more substituted amine than reactants
                            if checker.check_fg("Secondary amine", product) or checker.check_fg(
                                "Tertiary amine", product
                            ):
                                is_coupling = True
                                print("Found N-alkylation reaction based on functional groups")

                    # If we have a Boc group preserved through reaction, and it's a coupling reaction,
                    # then we have a protected amine coupling
                    if is_coupling:
                        found_pattern = True
                        print(f"Found protected amine coupling at depth {depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_pattern
