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
    This function detects a strategy involving esterification in the late stages of synthesis.
    In retrosynthesis, this includes both forward esterification and reverse hydrolysis reactions.
    """
    late_stage_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_esterification

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for esterification or hydrolysis reactions directly
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                ):
                    print(f"Detected esterification/hydrolysis reaction at depth {depth}")
                    late_stage_esterification = True
                    return

                # Fallback check using functional groups
                # Check for esterification (forward): carboxylic acid + alcohol → ester
                if checker.check_fg("Ester", product_smiles):
                    print(f"Product contains ester at depth {depth}")
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants_smiles
                    )

                    if has_carboxylic_acid and has_alcohol:
                        print(
                            f"Detected esterification pattern at depth {depth} (carboxylic acid + alcohol → ester)"
                        )
                        late_stage_esterification = True
                    elif has_carboxylic_acid:
                        print(
                            f"Detected possible esterification at depth {depth} (has carboxylic acid)"
                        )
                        late_stage_esterification = True

                # Check for hydrolysis (retrosynthetic): ester → carboxylic acid + alcohol
                if checker.check_fg("Carboxylic acid", product_smiles):
                    print(f"Product contains carboxylic acid at depth {depth}")
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)

                    if has_ester:
                        print(
                            f"Detected ester hydrolysis pattern at depth {depth} (ester → carboxylic acid)"
                        )
                        late_stage_esterification = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)
    print(f"Final result: late_stage_esterification = {late_stage_esterification}")

    return late_stage_esterification
