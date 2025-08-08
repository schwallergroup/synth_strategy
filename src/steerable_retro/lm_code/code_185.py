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
    This function detects if the synthesis includes a late-stage amidation
    (carboxylic acid to amide transformation in the last or second-to-last step).
    """
    has_late_amidation = False

    def dfs_traverse(node, depth=0, path=None):
        nonlocal has_late_amidation

        if path is None:
            path = []

        # Add current node to path
        path.append(node)

        # Check if this is a reaction node at depth 1 or 2 (late-stage)
        if node["type"] == "reaction" and 1 <= depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an amidation reaction using comprehensive list
                amidation_reactions = [
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Ester with primary amine to amide",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Schotten-Baumann to ester",
                    "Acyl chloride with primary amine to imide",
                ]

                # Check if any of the known amidation reactions match
                is_amidation = False
                for reaction in amidation_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        print(f"Matched reaction type: {reaction}")
                        is_amidation = True
                        break

                # Define functional groups to check
                acid_groups = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                amide_groups = ["Primary amide", "Secondary amide", "Tertiary amide"]

                # Check for acid groups in reactants
                has_acid = False
                for reactant in reactants:
                    for acid in acid_groups:
                        if checker.check_fg(acid, reactant):
                            print(f"Found {acid} in reactant: {reactant}")
                            has_acid = True

                # Check for amide groups in product
                has_amide = False
                for amide in amide_groups:
                    if checker.check_fg(amide, product):
                        print(f"Found {amide} in product: {product}")
                        has_amide = True

                # If we didn't match a known reaction but found acid → amide transformation
                if not is_amidation and has_acid and has_amide:
                    print("Detected amidation based on functional group transformation")
                    is_amidation = True

                if is_amidation and has_acid and has_amide:
                    has_late_amidation = True
                    print(f"Confirmed late-stage amidation at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            # Create a new path for each branch to avoid modifying the current path
            dfs_traverse(child, depth + 1, path.copy())

    dfs_traverse(route)
    print(
        f"Synthesis {'includes' if has_late_amidation else 'does not include'} late-stage amidation"
    )
    return has_late_amidation
