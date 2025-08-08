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
    This function detects a strategy where a morpholine substituent is preserved
    throughout the synthesis while other groups are modified.
    """
    # Track molecules with morpholine and their positions in the route
    molecules_with_morpholine = []
    # Track if we've found at least one reaction that modifies non-morpholine parts
    found_modifying_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_modifying_reaction

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if (
                mol and mol.GetNumAtoms() > 7
            ):  # Ensure molecule is large enough to have morpholine + other groups
                if checker.check_ring("morpholine", node["smiles"]):
                    molecules_with_morpholine.append((node["smiles"], depth))
                    print(f"Found morpholine in molecule: {node['smiles']} at depth {depth}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this reaction preserves morpholine while modifying other parts
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if both reactants and product contain morpholine
            reactants_with_morpholine = [
                r for r in reactants if checker.check_ring("morpholine", r)
            ]
            product_has_morpholine = checker.check_ring("morpholine", product)

            if reactants_with_morpholine and product_has_morpholine:
                # Verify that something else changed in the reaction
                if any(
                    Chem.MolFromSmiles(r).GetNumAtoms() != Chem.MolFromSmiles(product).GetNumAtoms()
                    for r in reactants_with_morpholine
                ):
                    found_modifying_reaction = True
                    print(
                        f"Found reaction that preserves morpholine while modifying other parts: {rsmi}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is detected if:
    # 1. We found multiple molecules with morpholine across different depths
    # 2. We found at least one reaction that modifies non-morpholine parts
    # 3. The target molecule (depth 0) has morpholine

    has_morpholine_at_different_depths = (
        len(set(depth for _, depth in molecules_with_morpholine)) > 1
    )
    target_has_morpholine = any(depth == 0 for _, depth in molecules_with_morpholine)

    strategy_detected = (
        len(molecules_with_morpholine) >= 2
        and has_morpholine_at_different_depths
        and target_has_morpholine
        and found_modifying_reaction
    )

    if strategy_detected:
        print("Detected morpholine-preserving strategy")
    else:
        if not molecules_with_morpholine:
            print("No morpholine found in any molecule")
        elif not has_morpholine_at_different_depths:
            print("Morpholine only found at a single depth")
        elif not target_has_morpholine:
            print("Target molecule does not contain morpholine")
        elif not found_modifying_reaction:
            print("No reactions found that preserve morpholine while modifying other parts")

    return strategy_detected
