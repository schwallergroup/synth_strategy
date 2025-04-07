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
    Detects if an amide formation occurs in the late stage of the synthesis.
    """
    late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation

        if node["type"] == "reaction" and depth <= 3:  # Late stage (depth ≤ 3)
            print(f"Examining reaction at depth {depth}")
            # Check if this is an amide formation
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # Check for amide formation reactions using the checker
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acyl chloride with ammonia to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Aminolysis of esters",
                    "Carboxylic acid to amide conversion",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide formation reaction: {reaction_type}")
                        late_amide_formation = True
                        break

                # If no specific reaction type matched, check for reactants and products
                if not late_amide_formation:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        print(f"Reactants: {reactants}")
                        print(f"Product: {product}")

                        # Check for amine in reactants
                        amine_reactants = []
                        for r in reactants:
                            if checker.check_fg("Primary amine", r) or checker.check_fg(
                                "Secondary amine", r
                            ):
                                amine_reactants.append(r)
                                print(f"Found amine in reactant: {r}")

                        # Check for acyl source in reactants
                        acyl_reactants = []
                        for r in reactants:
                            if (
                                checker.check_fg("Acyl halide", r)
                                or checker.check_fg("Carboxylic acid", r)
                                or checker.check_fg("Ester", r)
                                or checker.check_fg("Anhydride", r)
                            ):
                                acyl_reactants.append(r)
                                print(f"Found acyl source in reactant: {r}")

                        # Check for amide in product
                        amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                        product_amides = []
                        for amide_type in amide_types:
                            if checker.check_fg(amide_type, product):
                                product_amides.append(amide_type)
                                print(f"Found {amide_type} in product")

                        # Count amides in reactants
                        reactant_amides = []
                        for r in reactants:
                            for amide_type in amide_types:
                                if checker.check_fg(amide_type, r):
                                    reactant_amides.append((r, amide_type))
                                    print(f"Found {amide_type} in reactant: {r}")

                        # Verify that a new amide is formed
                        if (
                            len(amine_reactants) > 0
                            and len(acyl_reactants) > 0
                            and len(product_amides) > 0
                            and len(reactant_amides) < len(product_amides)
                        ):
                            print("Found amide formation based on functional group analysis")
                            late_amide_formation = True
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide formation: {late_amide_formation}")
    return late_amide_formation
