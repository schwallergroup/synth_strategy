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
    This function detects a synthetic strategy featuring mid-stage heterocycle formation
    followed by functional group interconversions with late-stage amide coupling.
    """
    # Initialize tracking variables
    heterocycle_formation_depth = None
    ester_hydrolysis_depth = None
    amide_coupling_depth = None
    reaction_sequence = []

    # List of common heterocycles to check
    heterocycles = [
        "imidazole",
        "pyrrole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "benzimidazole",
        "indole",
        "benzoxazole",
        "benzothiazole",
        "furan",
        "thiophene",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, ester_hydrolysis_depth, amide_coupling_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Store reaction information
            reaction_info = {
                "depth": depth,
                "rsmi": rsmi,
                "product_smiles": product_smiles,
                "reactants_smiles": reactants_smiles,
            }
            reaction_sequence.append(reaction_info)

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for heterocycle formation
            heterocycle_in_product = False
            heterocycle_in_reactants = False
            detected_heterocycle = None

            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, product_smiles):
                    heterocycle_in_product = True
                    detected_heterocycle = heterocycle
                    print(f"Found {heterocycle} in product at depth {depth}")

                    # Check if heterocycle was already in reactants
                    for reactant in reactants_smiles:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycle_in_reactants = True
                            print(f"Found {heterocycle} in reactant at depth {depth}")
                            break

                    if not heterocycle_in_reactants:
                        print(f"Detected {heterocycle} formation at depth {depth}")
                        heterocycle_formation_depth = depth
                        break

            # Check for ester hydrolysis or similar carboxylic acid formation
            if (
                checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                or checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )
                or checker.check_reaction("COOH ethyl deprotection", rsmi)
            ):
                print(f"Detected ester hydrolysis at depth {depth}")
                ester_hydrolysis_depth = depth
            # Additional check for ester to carboxylic acid conversion
            elif any(
                checker.check_fg("Ester", reactant) for reactant in reactants_smiles
            ) and checker.check_fg("Carboxylic acid", product_smiles):
                print(f"Detected ester to carboxylic acid conversion at depth {depth}")
                ester_hydrolysis_depth = depth

            # Check for amide coupling
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                or checker.check_reaction("{Schotten-Baumann_amide}", rsmi)
            ):
                print(f"Detected amide coupling at depth {depth}")
                amide_coupling_depth = depth
            # Additional check for amide formation
            elif not any(
                checker.check_fg("Primary amide", reactant)
                or checker.check_fg("Secondary amide", reactant)
                for reactant in reactants_smiles
            ) and (
                checker.check_fg("Primary amide", product_smiles)
                or checker.check_fg("Secondary amide", product_smiles)
            ):
                print(f"Detected amide formation at depth {depth}")
                amide_coupling_depth = depth

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required reactions were found
    all_reactions_found = (
        heterocycle_formation_depth is not None
        and ester_hydrolysis_depth is not None
        and amide_coupling_depth is not None
    )

    # Check if the reactions are in the correct order
    correct_order = False
    if all_reactions_found:
        # Amide coupling should be late-stage (lowest depth)
        # Ester hydrolysis and heterocycle formation can be in either order
        # as long as amide coupling is the latest stage
        correct_order = (
            amide_coupling_depth < ester_hydrolysis_depth
            and amide_coupling_depth < heterocycle_formation_depth
        )

        print(
            f"Reaction depths - Amide coupling: {amide_coupling_depth}, "
            f"Ester hydrolysis: {ester_hydrolysis_depth}, "
            f"Heterocycle formation: {heterocycle_formation_depth}"
        )
        print(f"Correct order: {correct_order}")

    # Strategy is present if all reactions are found and in the correct order
    strategy_present = all_reactions_found and correct_order

    if strategy_present:
        print("Detected heterocycle-centered linear synthesis with late-stage amide coupling")
    else:
        if not all_reactions_found:
            print("Not all required reactions were found")
        elif not correct_order:
            print("Reactions were not in the correct order")

    return strategy_present
