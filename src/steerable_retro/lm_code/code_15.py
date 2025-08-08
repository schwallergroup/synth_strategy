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
    This function detects a synthetic strategy involving late-stage amide bond formation
    in the final step of a linear synthesis.
    """
    amide_formation_at_depth_zero = False

    def dfs_traverse(node):
        nonlocal amide_formation_at_depth_zero

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract depth from ID field, default to 0 if not found (assume it's the final step)
            depth_match = None
            if "ID" in node["metadata"]:
                import re

                depth_match = re.search(r"Depth: (\d+)", node["metadata"]["ID"])
            depth = int(depth_match.group(1)) if depth_match else 0

            # Check if this is the final step (depth 0)
            if depth == 0:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Analyzing final step reaction: {rsmi}")
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if product contains an amide group
                    has_primary_amide = checker.check_fg("Primary amide", product_smiles)
                    has_secondary_amide = checker.check_fg("Secondary amide", product_smiles)
                    has_tertiary_amide = checker.check_fg("Tertiary amide", product_smiles)

                    if has_primary_amide:
                        print("Product contains primary amide")
                    if has_secondary_amide:
                        print("Product contains secondary amide")
                    if has_tertiary_amide:
                        print("Product contains tertiary amide")

                    if has_primary_amide or has_secondary_amide or has_tertiary_amide:
                        # Check if amide was formed in this reaction (not present in reactants)
                        amide_in_reactants = False
                        for reactant in reactants_smiles:
                            if (
                                checker.check_fg("Primary amide", reactant)
                                or checker.check_fg("Secondary amide", reactant)
                                or checker.check_fg("Tertiary amide", reactant)
                            ):
                                amide_in_reactants = True
                                print(f"Amide already present in reactant: {reactant}")
                                break

                        if not amide_in_reactants:
                            print("Amide not present in reactants - checking reaction types")
                            # Check for known amide formation reactions
                            amide_formation_reactions = [
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                "Carboxylic acid with primary amine to amide",
                                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                "Acyl chloride with secondary amine to amide",
                                "Acyl chloride with ammonia to amide",
                                "Ester with primary amine to amide",
                                "Ester with secondary amine to amide",
                                "Ester with ammonia to amide",
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                                "Schotten-Baumann to ester",
                                "{Schotten-Baumann_amide}",
                                "Acylation of primary amines",
                                "Acylation of secondary amines",
                                "Carboxylic acid to amide conversion",
                                "Nitrile and hydrogen peroxide to amide",
                            ]

                            for reaction_type in amide_formation_reactions:
                                if checker.check_reaction(reaction_type, rsmi):
                                    print(
                                        f"Detected amide formation in final step via {reaction_type}"
                                    )
                                    amide_formation_at_depth_zero = True
                                    break

                            # If no specific reaction detected, check for reactant functional groups
                            if not amide_formation_at_depth_zero:
                                print(
                                    "No specific amide formation reaction detected - checking reactant functional groups"
                                )
                                has_acid = False
                                has_amine = False
                                has_acyl_halide = False
                                has_ester = False
                                has_nitrile = False
                                has_anhydride = False

                                for reactant in reactants_smiles:
                                    if checker.check_fg("Carboxylic acid", reactant):
                                        has_acid = True
                                        print(f"Found carboxylic acid in reactant: {reactant}")
                                    if checker.check_fg(
                                        "Primary amine", reactant
                                    ) or checker.check_fg("Secondary amine", reactant):
                                        has_amine = True
                                        print(f"Found amine in reactant: {reactant}")
                                    if checker.check_fg("Acyl halide", reactant):
                                        has_acyl_halide = True
                                        print(f"Found acyl halide in reactant: {reactant}")
                                    if checker.check_fg("Ester", reactant):
                                        has_ester = True
                                        print(f"Found ester in reactant: {reactant}")
                                    if checker.check_fg("Nitrile", reactant):
                                        has_nitrile = True
                                        print(f"Found nitrile in reactant: {reactant}")
                                    if checker.check_fg("Anhydride", reactant):
                                        has_anhydride = True
                                        print(f"Found anhydride in reactant: {reactant}")

                                if (
                                    (has_acid and has_amine)
                                    or (has_acyl_halide and has_amine)
                                    or (has_ester and has_amine)
                                    or (has_nitrile)
                                    or (has_anhydride and has_amine)
                                ):
                                    print(
                                        "Detected amide formation in final step based on reactant functional groups"
                                    )
                                    amide_formation_at_depth_zero = True
                except Exception as e:
                    print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final result: {amide_formation_at_depth_zero}")
    return amide_formation_at_depth_zero
