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
    Detects if the synthesis route maintains a trifluoromethyl (CF3) group
    throughout the synthesis.
    """
    # Track if the target molecule has a CF3 group
    target_has_cf3 = False
    # Track if CF3 is preserved in all reactions
    cf3_preserved = True
    # Track if any starting material has CF3
    starting_materials_with_cf3 = False

    # Track all leaf nodes (potential starting materials)
    leaf_nodes = []

    def dfs_traverse(node, depth=0):
        nonlocal target_has_cf3, cf3_preserved, starting_materials_with_cf3

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)

            # Check if this is the target molecule (depth 0)
            if depth == 0:
                if has_cf3:
                    print(f"CF3 group found in target molecule: {mol_smiles}")
                    target_has_cf3 = True
                else:
                    print(f"Target molecule does not have CF3 group: {mol_smiles}")

            # Check if this is a starting material (leaf node or explicitly marked as in_stock)
            is_starting_material = node.get("in_stock", False) or len(node.get("children", [])) == 0

            if is_starting_material:
                leaf_nodes.append((mol_smiles, has_cf3, depth))
                if has_cf3:
                    print(f"CF3 group found in starting material at depth {depth}: {mol_smiles}")
                    starting_materials_with_cf3 = True
                else:
                    print(f"Starting material without CF3 at depth {depth}: {mol_smiles}")

            # Process reaction children (in retrosynthetic direction)
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    # Check if CF3 is preserved in this reaction
                    try:
                        rsmi = child["metadata"]["rsmi"]
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if product has CF3
                        product_has_cf3 = checker.check_fg("Trifluoro group", product)

                        # If product has CF3, at least one reactant should have CF3
                        if product_has_cf3:
                            reactants_have_cf3 = any(
                                checker.check_fg("Trifluoro group", r) for r in reactants
                            )

                            if not reactants_have_cf3:
                                print(
                                    f"CF3 group introduced in reaction at depth {depth+1}: {rsmi}"
                                )
                                cf3_preserved = False

                        # If product doesn't have CF3 but any reactant does, CF3 was lost
                        elif any(checker.check_fg("Trifluoro group", r) for r in reactants):
                            print(f"CF3 group lost in reaction at depth {depth+1}: {rsmi}")
                            cf3_preserved = False

                    except (KeyError, IndexError) as e:
                        print(f"Could not extract reaction information at depth {depth+1}: {e}")

                    # Continue traversal
                    dfs_traverse(child, depth + 1)

        elif node["type"] == "reaction":
            # Process molecule children (reactants in retrosynthetic direction)
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary of all leaf nodes found
    print("\nAll leaf nodes (potential starting materials):")
    for smiles, has_cf3, depth in leaf_nodes:
        status = "with CF3" if has_cf3 else "without CF3"
        print(f"Depth {depth}: {status} - {smiles}")

    # For a persistent CF3 strategy:
    # 1. Target molecule must have CF3
    # 2. CF3 must be preserved in all reactions (not introduced or removed)
    # 3. At least one starting material must have CF3

    # If we have CF3 in the target and it's preserved throughout, but no starting material
    # has CF3, then we need to check if there are any leaf nodes that might not be explicitly
    # marked as in_stock but are effectively starting materials
    if target_has_cf3 and cf3_preserved and not starting_materials_with_cf3:
        # Check if any leaf node has CF3
        for _, has_cf3, _ in leaf_nodes:
            if has_cf3:
                print("Found a leaf node with CF3 that might be a starting material")
                starting_materials_with_cf3 = True
                break

    result = target_has_cf3 and cf3_preserved and starting_materials_with_cf3
    print(f"\nTarget has CF3: {target_has_cf3}")
    print(f"CF3 preserved throughout: {cf3_preserved}")
    print(f"Starting materials with CF3: {starting_materials_with_cf3}")

    return result
