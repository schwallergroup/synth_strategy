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
    Detects if the synthesis maintains a methoxy-substituted aromatic ring
    throughout the entire synthetic route.
    """
    # Track continuous paths with methoxy groups
    continuous_paths = []

    def dfs_traverse(node, path=None, depth=0):
        if path is None:
            path = []

        current_path = path.copy()

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_methoxy = checker.check_fg("Ether", mol_smiles) and checker.check_ring(
                "benzene", mol_smiles
            )

            # For the target molecule (depth 0), initialize paths if methoxy is present
            if depth == 0 and has_methoxy:
                print(f"Found methoxy aromatic in target molecule: {mol_smiles}")
                current_path.append(True)
                continuous_paths.append(current_path)

            # For intermediate molecules, continue the path if methoxy is present
            elif depth > 0 and has_methoxy and len(path) > 0:
                print(f"Found methoxy aromatic at depth {depth}: {mol_smiles}")
                current_path.append(True)

            # If methoxy is not present, break the path
            elif depth > 0:
                print(f"No methoxy aromatic at depth {depth}: {mol_smiles}")
                current_path.append(False)

            # If this is a leaf node (starting material), check if the path is continuous
            if node.get("in_stock", False) and all(current_path):
                print(f"Found continuous methoxy path to starting material: {mol_smiles}")
                return True

        # For reaction nodes, just pass through the path
        elif node["type"] == "reaction":
            # Check if the reaction preserves the methoxy group
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")
                # We don't break the path for reaction nodes
                current_path.append(True)

        # Continue traversing children
        if len(node.get("children", [])) == 0 and all(current_path) and len(current_path) > 0:
            # Reached a leaf with a continuous methoxy path
            continuous_paths.append(current_path)
            return True

        result = False
        for child in node.get("children", []):
            if dfs_traverse(child, current_path, depth + 1):
                result = True

        return result

    # Start traversal
    has_continuous_path = dfs_traverse(route)

    # Check if we found any continuous paths with methoxy
    print(f"Continuous paths found: {len(continuous_paths)}")
    for i, path in enumerate(continuous_paths):
        print(f"Path {i+1}: {path}")

    return has_continuous_path or (
        len(continuous_paths) > 0 and any(all(path) for path in continuous_paths)
    )
