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
    Detects if a nitrogen heterocycle scaffold is preserved throughout the synthesis.
    """
    # List of nitrogen heterocycles to check
    n_heterocycles = [
        "pyridine",
        "pyrrole",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "piperidine",
        "pyrrolidine",
        "piperazine",
        "morpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "aziridine",
        "azetidine",
        "azepane",
    ]

    # Track which nitrogen heterocycles are present at each molecule node
    heterocycles_by_mol = {}
    mol_nodes_in_order = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            present_heterocycles = []

            # Check for each nitrogen heterocycle
            for ring_name in n_heterocycles:
                if checker.check_ring(ring_name, mol_smiles):
                    present_heterocycles.append(ring_name)
                    print(f"{ring_name} found in molecule: {mol_smiles}")

            # Store the heterocycles present in this molecule
            if present_heterocycles:
                heterocycles_by_mol[mol_smiles] = present_heterocycles
                mol_nodes_in_order.append(mol_smiles)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if at least one nitrogen heterocycle is found
    if not heterocycles_by_mol:
        print("No nitrogen heterocycles found in the synthesis route")
        return False

    # If we only have one molecule with heterocycles, can't determine preservation
    if len(mol_nodes_in_order) <= 1:
        print("Not enough molecules with nitrogen heterocycles to determine preservation")
        return False

    # Check for common heterocycles across all molecule nodes
    preserved_heterocycles = set(heterocycles_by_mol[mol_nodes_in_order[0]])

    for i in range(1, len(mol_nodes_in_order)):
        mol_smiles = mol_nodes_in_order[i]
        current_heterocycles = set(heterocycles_by_mol[mol_smiles])

        # Intersect with heterocycles in current molecule
        common_heterocycles = preserved_heterocycles & current_heterocycles

        # If no common heterocycles remain, preservation is broken
        if not common_heterocycles:
            print(
                f"No common nitrogen heterocycles preserved between {mol_nodes_in_order[i-1]} and {mol_smiles}"
            )
            return False

        preserved_heterocycles = common_heterocycles

    if preserved_heterocycles:
        print(
            f"Nitrogen heterocycles preserved throughout synthesis: {', '.join(preserved_heterocycles)}"
        )
        return True

    return False
