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
    This function detects if the synthetic route uses alcohol activation via mesylation,
    tosylation, triflation or similar strategies.
    """
    alcohol_activation = False

    def dfs_traverse(node):
        nonlocal alcohol_activation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            products = rsmi.split(">")[-1]

            # Check for alcohol in reactants
            if (
                checker.check_fg("Primary alcohol", reactants)
                or checker.check_fg("Secondary alcohol", reactants)
                or checker.check_fg("Tertiary alcohol", reactants)
                or checker.check_fg("Aromatic alcohol", reactants)
            ):

                # Check for mesylation (Formation of Sulfonic Esters)
                if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                    print(f"Found alcohol activation via mesylation: {rsmi}")
                    alcohol_activation = True

                # Check for tosylation
                elif checker.check_fg("Tosylate", products):
                    print(f"Found alcohol activation via tosylation: {rsmi}")
                    alcohol_activation = True

                # Check for triflation
                elif checker.check_fg("Triflate", products):
                    print(f"Found alcohol activation via triflation: {rsmi}")
                    alcohol_activation = True

                # Check for other sulfonates
                elif checker.check_fg("Mesylate", products):
                    print(f"Found alcohol activation via mesylation: {rsmi}")
                    alcohol_activation = True

                # Check for conversion to alkyl halides (another activation strategy)
                elif (
                    checker.check_fg("Primary halide", products)
                    or checker.check_fg("Secondary halide", products)
                    or checker.check_fg("Tertiary halide", products)
                ) and (
                    checker.check_reaction("Alkyl chlorides from alcohols", rsmi)
                    or checker.check_reaction("Alkyl bromides from alcohols", rsmi)
                    or checker.check_reaction("Alkyl iodides from alcohols", rsmi)
                ):
                    print(f"Found alcohol activation via halogenation: {rsmi}")
                    alcohol_activation = True

                # Check for alcohol to triflate conversion
                elif checker.check_reaction("Alcohol to triflate conversion", rsmi):
                    print(f"Found alcohol activation via triflation: {rsmi}")
                    alcohol_activation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return alcohol_activation
