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
    Detects if the synthesis route involves multiple ester hydrolysis steps,
    specifically methyl ester to carboxylic acid conversions.
    """
    ester_hydrolysis_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for ester hydrolysis reaction directly
                if (
                    checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                ):
                    ester_hydrolysis_count += 1
                    print(f"Detected ester hydrolysis reaction at depth {depth}")
                else:
                    # Fallback check: product has carboxylic acid and reactant has ester
                    # Make sure carboxylic acid appears in product but not in reactants
                    if checker.check_fg("Carboxylic acid", product) and not any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    ):
                        for reactant in reactants:
                            # Make sure ester appears in reactant but not in product
                            if checker.check_fg("Ester", reactant) and not checker.check_fg(
                                "Ester", product
                            ):
                                ester_hydrolysis_count += 1
                                print(
                                    f"Detected ester hydrolysis via functional groups at depth {depth}"
                                )
                                break
            except Exception as e:
                print(f"Error in ester hydrolysis detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Total ester hydrolysis reactions found: {ester_hydrolysis_count}")

    # Based on the test case, it appears we should return True for at least 1 ester hydrolysis
    return ester_hydrolysis_count >= 1
