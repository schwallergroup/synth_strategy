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
    Detects synthesis routes that involve multiple nitrogen heterocycles.
    """
    heterocycle_count = 0
    detected_heterocycles = set()

    # List of nitrogen-containing heterocycles to check
    nitrogen_heterocycles = [
        "pyrazole",
        "pyridine",
        "triazole",
        "imidazole",
        "piperidine",
        "piperazine",
        "morpholine",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "pyrrolidine",
        "azetidine",
        "aziridine",
        "azepane",
        "diazepane",
        "oxazole",
        "thiazole",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_count, detected_heterocycles

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for each nitrogen heterocycle
            for heterocycle in nitrogen_heterocycles:
                if (
                    checker.check_ring(heterocycle, mol_smiles)
                    and heterocycle not in detected_heterocycles
                ):
                    detected_heterocycles.add(heterocycle)
                    heterocycle_count += 1
                    print(f"Found {heterocycle} in molecule: {mol_smiles}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = heterocycle_count >= 2
    print(
        f"Heterocycle-rich synthesis: {result} (count: {heterocycle_count}, detected: {detected_heterocycles})"
    )
    return result
