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
    Detects if heterocycles are preserved throughout the synthesis route.
    """

    def get_heterocycles(mol_smiles):
        heterocycles = []
        for ring_name in [
            "furan",
            "pyrrole",
            "pyridine",
            "pyrazole",
            "imidazole",
            "oxazole",
            "thiazole",
            "pyrimidine",
            "pyrazine",
            "pyridazine",
            "triazole",
            "tetrazole",
        ]:
            if checker.check_ring(ring_name, mol_smiles):
                heterocycles.append(ring_name)
        return heterocycles

    # First identify target heterocycles
    target_heterocycles = []
    if route["type"] == "mol":
        target_heterocycles = get_heterocycles(route["smiles"])
        print(f"Target molecule contains heterocycles: {target_heterocycles}")

    # If no heterocycles in target, strategy doesn't apply
    if not target_heterocycles:
        print("No heterocycles in target molecule")
        return True

    def dfs(node, depth=0):
        if node["type"] == "mol" and not node.get("in_stock", False):
            current_heterocycles = get_heterocycles(node["smiles"])

            # Check if all target heterocycles are present
            for cycle in target_heterocycles:
                if cycle not in current_heterocycles:
                    print(f"Heterocycle {cycle} not preserved at depth {depth}")
                    return False

        # Continue DFS traversal
        all_preserved = True
        for child in node.get("children", []):
            if not dfs(child, depth + 1):
                all_preserved = False

        return all_preserved

    return dfs(route)
