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
    This function detects if the synthetic route involves amide bond formation
    in the final step (depth 1).
    """
    late_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_coupling

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # For reaction nodes, check if it's an amide coupling
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check for amide formation at depth 1 (first reaction in retrosynthetic route)
            if depth == 1:
                print("Checking for late-stage amide coupling...")

                # Check if this is an amide coupling reaction using predefined reaction types
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                    "Carboxylic acid to amide conversion",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                ]

                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Late-stage amide coupling detected: {reaction_type}")
                        late_amide_coupling = True
                        return

                # If no specific reaction type matched, check for functional group changes
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reactants: {reactants}")
                    print(f"Analyzing product: {product}")

                    # Check for carboxylic acid derivatives and amine in reactants
                    carboxylic_acid_present = False
                    amine_present = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            carboxylic_acid_present = True
                            print(f"Found carboxylic acid derivative in reactant: {reactant}")

                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            amine_present = True
                            print(f"Found amine in reactant: {reactant}")

                    # Check if product has amide
                    product_has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if product_has_amide:
                        print(f"Found amide in product: {product}")

                    # Check if any reactant already has the amide
                    reactants_have_amide = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                        ):
                            reactants_have_amide = True
                            print(f"Found amide in reactant: {reactant}")
                            break

                    # Confirm a new amide is formed
                    if (
                        product_has_amide
                        and not reactants_have_amide
                        and carboxylic_acid_present
                        and amine_present
                    ):
                        print(
                            "Late-stage amide coupling detected through functional group analysis"
                        )
                        late_amide_coupling = True
                except Exception as e:
                    print(f"Error in functional group analysis: {e}")

        # Continue traversing
        for child in node.get("children", []):
            if not late_amide_coupling:  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {late_amide_coupling}")
    return late_amide_coupling
