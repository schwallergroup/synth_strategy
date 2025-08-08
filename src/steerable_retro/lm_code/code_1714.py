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
    This function detects if bromination occurs early in the synthesis and the bromine
    is maintained until late stages.
    """
    # Track depths where bromine is present and where bromination occurs
    depths_with_bromine = set()
    depths_with_bromination = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            # Check for bromine-containing functional groups
            if (
                checker.check_fg("Aromatic halide", node["smiles"])
                or checker.check_fg("Primary halide", node["smiles"])
                or checker.check_fg("Secondary halide", node["smiles"])
                or checker.check_fg("Tertiary halide", node["smiles"])
                or checker.check_fg("Alkenyl halide", node["smiles"])
            ):

                # Verify it's specifically bromine by checking SMARTS
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[Br]")):
                    depths_with_bromine.add(depth)
                    print(f"Found bromine at depth {depth}: {node['smiles']}")

        # Check for bromination reactions
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            if checker.check_reaction("Aromatic bromination", rxn_smiles) or checker.check_reaction(
                "Bromination", rxn_smiles
            ):
                depths_with_bromination.add(depth)
                print(f"Found bromination reaction at depth {depth}: {rxn_smiles}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Early stage is considered as depth > max_depth/2
    # Late stage is considered as depth <= 2
    early_stage_threshold = max_depth // 2

    # Check conditions:
    # 1. We have bromination reactions in early stages
    # 2. Bromine is present in both early and late stages
    has_early_bromination = any(d >= early_stage_threshold for d in depths_with_bromination)
    has_early_bromine = any(d >= early_stage_threshold for d in depths_with_bromine)
    has_late_bromine = any(d <= 2 for d in depths_with_bromine)

    print(f"Max depth: {max_depth}")
    print(f"Depths with bromine: {sorted(depths_with_bromine)}")
    print(f"Depths with bromination: {sorted(depths_with_bromination)}")
    print(f"Early stage threshold: {early_stage_threshold}")
    print(f"Has early bromination: {has_early_bromination}")
    print(f"Has early bromine: {has_early_bromine}")
    print(f"Has late bromine: {has_late_bromine}")

    # Return True if we have bromination in early stages and bromine is maintained until late stages
    if max_depth > 2 and (has_early_bromination or has_early_bromine) and has_late_bromine:
        print("Bromine is introduced early and maintained until late stages")
        return True

    return False
