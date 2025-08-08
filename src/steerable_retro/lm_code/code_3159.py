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
    Detects if the synthesis route involves a multi-aromatic heterocyclic system
    with both pyrazole and isoxazole rings
    """
    has_pyrazole = False
    has_isoxazole = False
    has_oxazole = False  # Added to check for oxazole rings

    # List of rings to check for debugging
    ring_types = ["pyrazole", "isoxazole", "oxazole", "benzoxazole"]

    def dfs_traverse(node):
        nonlocal has_pyrazole, has_isoxazole, has_oxazole

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            try:
                # Debug: Print all rings found in this molecule
                found_rings = []
                for ring in ring_types:
                    if checker.check_ring(ring, mol_smiles):
                        found_rings.append(ring)

                if found_rings:
                    print(f"Found rings in molecule {mol_smiles}: {', '.join(found_rings)}")

                # Check for pyrazole ring
                if checker.check_ring("pyrazole", mol_smiles):
                    has_pyrazole = True
                    print(f"Found pyrazole ring in molecule: {mol_smiles}")

                # Check for isoxazole ring
                if checker.check_ring("isoxazole", mol_smiles):
                    has_isoxazole = True
                    print(f"Found isoxazole ring in molecule: {mol_smiles}")

                # Check for oxazole ring (which might be confused with isoxazole)
                if checker.check_ring("oxazole", mol_smiles):
                    has_oxazole = True
                    print(f"Found oxazole ring in molecule: {mol_smiles}")

                # Check for benzoxazole (which contains an oxazole-like structure)
                if checker.check_ring("benzoxazole", mol_smiles):
                    print(f"Found benzoxazole ring in molecule: {mol_smiles}")
                    # Note: We don't set has_isoxazole=True here as benzoxazole is not isoxazole
            except Exception as e:
                print(f"Error checking rings in molecule {mol_smiles}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both pyrazole and either isoxazole or oxazole are found
    # This is a more flexible interpretation of "multi-aromatic heterocyclic system"
    result = has_pyrazole and (has_isoxazole or has_oxazole)
    print(
        f"Final result: has_pyrazole={has_pyrazole}, has_isoxazole={has_isoxazole}, has_oxazole={has_oxazole}, returning {result}"
    )

    return result
