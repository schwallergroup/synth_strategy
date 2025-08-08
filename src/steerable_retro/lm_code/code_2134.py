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
    This function detects late-stage amide formation in the synthetic route.
    Late-stage is defined as reactions occurring at depths 0-2 in the synthesis tree.
    """
    found_late_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_amide

        # Check reaction nodes at depths 0-2 (late-stage)
        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if this is an amide formation reaction using the checker
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
                ]

                # Check if any of the amide formation reactions match
                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Found late-stage amide formation reaction: {reaction_type} at depth {depth}"
                        )
                        found_late_amide = True
                        return

                # If no specific reaction matched, check for functional group transformation
                # Split reactants and check product
                reactants = reactants_str.split(".")
                product = product_str

                # Check if product contains amide
                if (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                ):

                    # Check if reactants contain amide precursors
                    has_acid_precursor = False
                    has_amine = False

                    for reactant in reactants:
                        # Check for carboxylic acid or derivatives
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            has_acid_precursor = True

                        # Check for amine components
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Ammonia", reactant)
                        ):
                            has_amine = True

                    # If we have both precursors and the product has an amide, it's likely amide formation
                    if has_acid_precursor and has_amine:
                        print(
                            f"Found late-stage amide formation (functional group analysis) at depth {depth}"
                        )
                        found_late_amide = True

        # Process children
        for child in node.get("children", []):
            if not found_late_amide:  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide formation strategy detected: {found_late_amide}")
    return found_late_amide
