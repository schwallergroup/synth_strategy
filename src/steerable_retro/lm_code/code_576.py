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
    Detects if the synthesis route involves amide formation in the final step.
    """
    result = False

    # Find the final synthetic step (the first reaction node in retrosynthetic analysis)
    final_step = None

    def find_final_step(node):
        nonlocal final_step
        if node["type"] == "mol" and not node.get("in_stock", False):
            # This is the target molecule
            if node.get("children", []):
                for child in node["children"]:
                    if child["type"] == "reaction":
                        final_step = child
                        return True
        return False

    # First identify the final step
    find_final_step(route)

    print(f"Final step identified: {final_step is not None}")

    if final_step:
        try:
            # Check if this is an amide formation reaction using reaction checker
            if "rsmi" in final_step.get("metadata", {}):
                rsmi = final_step["metadata"]["rsmi"]
                print(f"Final step reaction SMILES: {rsmi}")

                # Check for known amide formation reaction types
                amide_reaction_types = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Schotten-Baumann to ester",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                for rxn_type in amide_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected amide formation reaction: {rxn_type}")
                        result = True
                        break

                # If no specific reaction type matched, check for functional group transformation
                if not result:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for amine and carboxylic acid in reactants
                    has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants)

                    # Check for amide in product
                    has_primary_amide = checker.check_fg("Primary amide", product)
                    has_secondary_amide = checker.check_fg("Secondary amide", product)
                    has_tertiary_amide = checker.check_fg("Tertiary amide", product)

                    print(
                        f"Reactants - Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}"
                    )
                    print(
                        f"Reactants - Carboxylic acid: {has_carboxylic_acid}, Acyl halide: {has_acyl_halide}, Ester: {has_ester}, Anhydride: {has_anhydride}"
                    )
                    print(
                        f"Product - Primary amide: {has_primary_amide}, Secondary amide: {has_secondary_amide}, Tertiary amide: {has_tertiary_amide}"
                    )

                    # Check for amide formation conditions
                    if (
                        (has_primary_amine or has_secondary_amine)
                        and (has_carboxylic_acid or has_acyl_halide or has_ester or has_anhydride)
                        and (has_primary_amide or has_secondary_amide or has_tertiary_amide)
                    ):
                        print(
                            "Detected amide formation in final step based on functional group analysis"
                        )
                        result = True
        except Exception as e:
            print(f"Error analyzing reaction: {e}")

    return result
