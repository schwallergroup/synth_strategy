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
    Detects if the synthesis route involves a convergent strategy with a piperazine scaffold.
    This means multiple branches converge to form a molecule containing piperazine.
    """
    # Track if we found piperazine
    found_piperazine = False
    # Track if we found a convergent reaction (multiple reactants)
    found_convergent = False

    def dfs(node, depth=0):
        nonlocal found_piperazine, found_convergent

        if node["type"] == "mol":
            # Check if this molecule contains piperazine
            if checker.check_ring("piperazine", node["smiles"]):
                found_piperazine = True
                print(f"Found piperazine in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                # Split into reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a convergent reaction (multiple reactants)
                if len(reactants) > 1:
                    found_convergent = True
                    print(f"Found convergent reaction with {len(reactants)} reactants: {rsmi}")

                # Check if product contains piperazine
                if checker.check_ring("piperazine", product):
                    found_piperazine = True
                    print(f"Found piperazine in product: {product}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    result = found_piperazine and found_convergent
    print(f"Convergent synthesis with piperazine: {result}")
    return result
