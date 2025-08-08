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
    This function detects if the route uses a late-stage amide formation strategy
    (amide formation occurs in the last 2 steps of the synthesis).
    """
    has_late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_formation

        # Process current node
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is an amide formation reaction using the checker function
            amide_formation_reactions = [
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
                "Carboxylic acid to amide conversion",
                "Acylation of primary amines",
                "Acylation of secondary amines",
            ]

            is_amide_formation = any(
                checker.check_reaction(rxn_type, rsmi) for rxn_type in amide_formation_reactions
            )

            # If not detected by reaction type, check for amide formation by examining product and reactants
            if not is_amide_formation:
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if amide is in product
                has_amide_in_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                # Check if amide is not in any reactant
                has_amide_in_reactants = False
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                    ):
                        has_amide_in_reactants = True
                        break

                # Check for reactants that could form amides
                has_acid_or_acyl = False
                has_amine = False

                for reactant in reactants:
                    if (
                        checker.check_fg("Carboxylic acid", reactant)
                        or checker.check_fg("Acyl halide", reactant)
                        or checker.check_fg("Ester", reactant)
                        or checker.check_fg("Anhydride", reactant)
                    ):
                        has_acid_or_acyl = True

                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        has_amine = True

                # Amide formation if: amide in product + not in reactants + appropriate reactants present
                is_amide_formation = (
                    has_amide_in_product
                    and not has_amide_in_reactants
                    and has_acid_or_acyl
                    and has_amine
                )

                if is_amide_formation:
                    print(f"Detected amide formation by functional group analysis at depth {depth}")
                    print(f"  Product has amide: {has_amide_in_product}")
                    print(f"  Reactants have amide: {has_amide_in_reactants}")
                    print(f"  Reactants have acid/acyl: {has_acid_or_acyl}")
                    print(f"  Reactants have amine: {has_amine}")

            # Check if this is a late-stage reaction (depth <= 1 for last 2 steps)
            if is_amide_formation and depth <= 1:
                has_late_amide_formation = True
                print(f"Late-stage amide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # For reaction nodes, increment depth when moving to children
            next_depth = depth
            if node["type"] == "reaction":
                next_depth = depth + 1
            dfs_traverse(child, next_depth)

    # Start traversal from depth 0 (target molecule)
    dfs_traverse(route)
    print(f"Late-stage amide formation: {has_late_amide_formation}")
    return has_late_amide_formation
