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
    This function detects a synthetic strategy involving pyrazole ring formation
    from nitrile-containing precursors and hydrazine or related compounds.
    """
    pyrazole_formed = False

    def dfs_traverse(node):
        nonlocal pyrazole_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for pyrazole in product
            pyrazole_in_product = checker.check_ring("pyrazole", product_smiles)

            # Check if pyrazole was not in reactants
            pyrazole_in_reactants = any(checker.check_ring("pyrazole", r) for r in reactants_smiles)

            # Only proceed if we're forming a pyrazole (not present in reactants but present in product)
            if pyrazole_in_product and not pyrazole_in_reactants:
                print(
                    f"Potential pyrazole formation detected - pyrazole in product but not in reactants"
                )

                # Check for various pyrazole formation reactions
                is_pyrazole_reaction = (
                    checker.check_reaction("pyrazole", rsmi)
                    or checker.check_reaction("{pyrazole}", rsmi)
                    or checker.check_reaction("Pyrazole formation", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkyne", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkene", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkyne", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkene", rsmi)
                )

                # Check for nitrile precursors in reactants
                nitrile_in_reactants = any(checker.check_fg("Nitrile", r) for r in reactants_smiles)

                # Check for hydrazine or related compounds in reactants
                hydrazine_related_in_reactants = any(
                    checker.check_fg("Hydrazine", r)
                    or checker.check_fg("Acylhydrazine", r)
                    or checker.check_fg("Hydrazone", r)
                    or checker.check_fg("Hydrazone amide", r)
                    or checker.check_fg("Diazene", r)
                    or checker.check_fg("Diazo", r)
                    for r in reactants_smiles
                )

                print(f"Is pyrazole reaction: {is_pyrazole_reaction}")
                print(f"Nitrile in reactants: {nitrile_in_reactants}")
                print(f"Hydrazine or related in reactants: {hydrazine_related_in_reactants}")

                # Check if this is a pyrazole formation from nitrile strategy
                if is_pyrazole_reaction or (
                    nitrile_in_reactants and hydrazine_related_in_reactants
                ):
                    print("Pyrazole ring formation from nitrile strategy detected!")
                    pyrazole_formed = True

            # Also check for other reactions that might be part of a multi-step pyrazole formation
            elif not pyrazole_formed:
                # Check for hydrazone formation (intermediate step)
                if checker.check_reaction("Ketone/aldehyde to hydrazone", rsmi):
                    print(
                        "Hydrazone formation detected - potential intermediate for pyrazole synthesis"
                    )

                # Check for diazoalkane formation (intermediate step)
                if checker.check_reaction("Hydrazone oxidation to diazoalkane", rsmi):
                    print(
                        "Diazoalkane formation detected - potential intermediate for pyrazole synthesis"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Pyrazole formation from nitrile strategy detected: {pyrazole_formed}")
    return pyrazole_formed
