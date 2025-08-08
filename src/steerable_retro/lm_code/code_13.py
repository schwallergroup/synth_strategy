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
    This function detects a synthetic strategy involving transformations on
    a substrate containing both halogenated aromatic rings and amide groups.
    """
    halogenated_amide_present = False
    multiple_transformations = 0

    def dfs_traverse(node):
        nonlocal halogenated_amide_present, multiple_transformations

        if node["type"] == "mol":
            # Check if molecule contains both halogenated aromatic and amide
            if "smiles" in node:
                mol_smiles = node["smiles"]
                has_halogen = checker.check_fg("Aromatic halide", mol_smiles)
                has_amide = (
                    checker.check_fg("Primary amide", mol_smiles)
                    or checker.check_fg("Secondary amide", mol_smiles)
                    or checker.check_fg("Tertiary amide", mol_smiles)
                )

                if has_halogen and has_amide:
                    halogenated_amide_present = True
                    print(f"Found molecule with halogenated aromatic and amide: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has both halogenated aromatic and amide
                product_has_halogen = checker.check_fg("Aromatic halide", product)
                product_has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                if product_has_halogen and product_has_amide:
                    # Check if any reactant also has both groups
                    for reactant in reactants:
                        reactant_has_halogen = checker.check_fg("Aromatic halide", reactant)
                        reactant_has_amide = (
                            checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                        )

                        if reactant_has_halogen and reactant_has_amide:
                            multiple_transformations += 1
                            print(
                                f"Detected transformation on halogenated amide substrate (count: {multiple_transformations})"
                            )
                            print(f"Reaction: {rsmi}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if halogenated amide is present and undergoes multiple transformations
    result = halogenated_amide_present and multiple_transformations >= 2
    print(
        f"Final result: halogenated_amide_present={halogenated_amide_present}, multiple_transformations={multiple_transformations}"
    )
    return result
