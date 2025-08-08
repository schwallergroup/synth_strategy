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
    Detects if the synthesis route involves the conversion of a hydroxyl group to a chloride.
    """
    conversion_detected = False

    def dfs_traverse(node):
        nonlocal conversion_detected

        if conversion_detected:
            return  # Early return if already detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for specific alcohol to chloride reactions
                alcohol_to_chloride_reactions = [
                    "Alcohol to chloride_sulfonyl chloride",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_CHCl3",
                    "Alcohol to chloride_CH2Cl2",
                    "Alcohol to chloride_PCl5_ortho",
                    "Alcohol to chloride_POCl3_ortho",
                    "Alcohol to chloride_POCl3_para",
                    "Alcohol to chloride_POCl3",
                    "Alcohol to chloride_HCl",
                    "Alcohol to chloride_Salt",
                    "Alcohol to chloride_Other",
                ]

                for reaction_type in alcohol_to_chloride_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected hydroxyl to chloride conversion: {reaction_type}")
                        conversion_detected = True
                        return

                # If specific reaction check failed, try checking functional group conversion
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chloride in product
                if (
                    checker.check_fg("Primary halide", product)
                    or checker.check_fg("Secondary halide", product)
                    or checker.check_fg("Tertiary halide", product)
                    or checker.check_fg("Aromatic halide", product)
                ):

                    # Check for hydroxyl in reactants
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                        ):

                            # This is a potential hydroxyl to chloride conversion
                            # Ideally, we would check atom mapping to confirm the exact position
                            print(f"Detected potential hydroxyl to chloride conversion")
                            conversion_detected = True
                            return

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return conversion_detected
