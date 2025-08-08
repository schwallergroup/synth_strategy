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
    This function detects a synthetic strategy where aromatic scaffolds
    (pyrazole and chlorobenzene) are preserved throughout the synthesis.
    """
    # Track presence of key scaffolds at each step
    scaffold_tracking = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node.get("smiles"):
            smiles = node["smiles"]

            # Check for pyrazole ring and chlorobenzene (aromatic halide)
            has_pyrazole = checker.check_ring("pyrazole", smiles)
            has_chlorobenzene = checker.check_fg("Aromatic halide", smiles)

            scaffold_tracking.append(
                {
                    "depth": depth,
                    "smiles": smiles,
                    "pyrazole": has_pyrazole,
                    "chlorobenzene": has_chlorobenzene,
                }
            )

            print(
                f"Depth {depth}: Molecule {smiles} - Pyrazole: {has_pyrazole}, Chlorobenzene: {has_chlorobenzene}"
            )

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}: Reaction {rsmi}")

            # Check if this reaction preserves the scaffolds
            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants and products contain the scaffolds
                reactants_pyrazole = any(
                    checker.check_ring("pyrazole", r) for r in reactants_smiles
                )
                reactants_chlorobenzene = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                )

                product_pyrazole = checker.check_ring("pyrazole", product_smiles)
                product_chlorobenzene = checker.check_fg("Aromatic halide", product_smiles)

                # If scaffolds are in reactants, they should also be in products for preservation
                if (reactants_pyrazole and not product_pyrazole) or (
                    reactants_chlorobenzene and not product_chlorobenzene
                ):
                    print(f"Scaffold not preserved in reaction: {rsmi}")
                    return False
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            result = dfs_traverse(child, depth + 1)
            if result is False:  # If any reaction doesn't preserve scaffolds, return False
                return False

    # Start traversal from root
    dfs_traverse(route)

    # Analyze scaffold tracking data
    if not scaffold_tracking:
        print("No molecules found in the route")
        return False

    # Check if both scaffolds are present in at least one molecule
    any_pyrazole = any(entry["pyrazole"] for entry in scaffold_tracking)
    any_chlorobenzene = any(entry["chlorobenzene"] for entry in scaffold_tracking)

    if not (any_pyrazole and any_chlorobenzene):
        print("One or both scaffolds not found in any molecule")
        return False

    # Group by depth to analyze each synthesis stage
    by_depth = {}
    for entry in scaffold_tracking:
        depth = entry["depth"]
        if depth not in by_depth:
            by_depth[depth] = []
        by_depth[depth].append(entry)

    # Check if both scaffolds are present at each depth where they should be
    # For scaffold preservation, if a scaffold appears at any point, it should persist
    max_depth = max(by_depth.keys()) if by_depth else 0
    pyrazole_appears_at = None
    chlorobenzene_appears_at = None

    for depth in range(max_depth + 1):
        if depth not in by_depth:
            continue

        depth_entries = by_depth[depth]
        depth_has_pyrazole = any(entry["pyrazole"] for entry in depth_entries)
        depth_has_chlorobenzene = any(entry["chlorobenzene"] for entry in depth_entries)

        # Track when scaffolds first appear
        if depth_has_pyrazole and pyrazole_appears_at is None:
            pyrazole_appears_at = depth
        if depth_has_chlorobenzene and chlorobenzene_appears_at is None:
            chlorobenzene_appears_at = depth

        # Once a scaffold appears, it should persist in later stages (lower depths)
        if (
            pyrazole_appears_at is not None
            and depth < pyrazole_appears_at
            and not depth_has_pyrazole
        ):
            print(f"Pyrazole not preserved at depth {depth}")
            return False
        if (
            chlorobenzene_appears_at is not None
            and depth < chlorobenzene_appears_at
            and not depth_has_chlorobenzene
        ):
            print(f"Chlorobenzene not preserved at depth {depth}")
            return False

    # Both scaffolds must appear and be preserved
    if pyrazole_appears_at is None or chlorobenzene_appears_at is None:
        print("One or both scaffolds never appear in the synthesis")
        return False

    print("Detected aromatic scaffold preservation strategy")
    return True
