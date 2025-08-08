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
    Detects if the synthetic route involves amide formation as the final step.
    """
    final_step_amide = False
    first_reaction_found = False

    def dfs_traverse(node):
        nonlocal final_step_amide, first_reaction_found

        # If we already found the first reaction, no need to process further nodes
        if first_reaction_found:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Mark this as the first reaction encountered (final step in retrosynthesis)
            first_reaction_found = True
            rsmi = node["metadata"]["rsmi"]

            print(f"Analyzing first reaction: {rsmi}")

            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Product: {product}")
                print(f"Reactants: {reactants}")

                # Check if this is an amide formation reaction using reaction checkers
                is_amide_formation = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                )

                print(f"Is amide formation reaction: {is_amide_formation}")

                # Check if product has amide
                product_has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                print(f"Product has amide: {product_has_amide}")

                # Check if reactants have necessary groups for amide formation
                has_amine = False
                has_acid_or_derivative = False
                has_ammonia = False

                for r in reactants:
                    # Check for amines
                    if checker.check_fg("Primary amine", r) or checker.check_fg(
                        "Secondary amine", r
                    ):
                        has_amine = True
                        print(f"Found amine in reactant: {r}")

                    # Check for ammonia specifically
                    if "[NH3" in r or r.strip() == "N":
                        has_ammonia = True
                        has_amine = True
                        print(f"Found ammonia in reactant: {r}")

                    # Check for acid derivatives
                    if (
                        checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Ester", r)
                        or checker.check_fg("Anhydride", r)
                        or "CO[C" in r
                    ):  # Additional check for methyl ester
                        has_acid_or_derivative = True
                        print(f"Found acid/derivative in reactant: {r}")

                # Check if amide is newly formed (not present in reactants)
                reactants_have_amide = False
                for r in reactants:
                    if (
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                    ):
                        reactants_have_amide = True
                        print(f"Found amide in reactant: {r}")

                # Check for ester to primary amide conversion (common with ammonia)
                ester_to_amide = False
                if has_ammonia and any("CO[C" in r for r in reactants) and product_has_amide:
                    ester_to_amide = True
                    print("Detected ester to primary amide conversion with ammonia")

                # Determine if this is an amide formation
                if is_amide_formation:
                    final_step_amide = True
                    print(f"Confirmed amide formation via reaction type check: {rsmi}")
                elif (
                    product_has_amide
                    and (has_amine or has_ammonia)
                    and has_acid_or_derivative
                    and not reactants_have_amide
                ):
                    final_step_amide = True
                    print(f"Confirmed amide formation via functional group analysis: {rsmi}")
                elif ester_to_amide:
                    final_step_amide = True
                    print(f"Confirmed ester to amide conversion: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction SMILES: {rsmi}, Error: {e}")

        # Process children only if we haven't found the first reaction yet
        if not first_reaction_found:
            for child in node.get("children", []):
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return final_step_amide
