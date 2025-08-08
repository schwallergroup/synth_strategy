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
    Detects synthesis routes that utilize methoxy-substituted aromatic building blocks,
    which are common in medicinal chemistry and can serve as handles for further functionalization.
    """
    methoxy_aromatic_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal methoxy_aromatic_count

        print(
            f"Traversing node at depth {depth}: {node.get('type', 'unknown')} - {node.get('smiles', 'no smiles')}"
        )

        if node.get("type") == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                mol = Chem.MolFromSmiles(mol_smiles)

                if mol:
                    # Check for methoxy group on aromatic ring
                    if checker.check_fg("Ether", mol_smiles) and mol.HasSubstructMatch(
                        Chem.MolFromSmarts("c[O][CH3]")
                    ):
                        methoxy_aromatic_count += 1
                        print(f"Found methoxy aromatic building block: {mol_smiles}")
            except Exception as e:
                print(f"Error processing molecule {node.get('smiles', 'unknown')}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Found {methoxy_aromatic_count} methoxy aromatic building blocks")

    # Return True if at least one methoxy aromatic building block is used
    return methoxy_aromatic_count > 0
