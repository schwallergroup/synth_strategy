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
    Detects convergent synthesis with thioether linkage formation and pyrimidinone construction.
    """
    # Initialize tracking variables
    has_thioether_formation = False
    has_pyrimidinone_formation = False
    has_convergent_pattern = False
    depth_2_fragments = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_thioether_formation, has_pyrimidinone_formation, has_convergent_pattern, depth_2_fragments

        # Count fragments at depth 2
        if depth == 2 and node["type"] == "mol" and not node.get("in_stock", False):
            depth_2_fragments += 1
            print(f"Found fragment at depth 2: {node['smiles']}")

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thioether formation (C-S-C)
            if depth <= 2:
                product_has_thioether = checker.check_fg("Monosulfide", product)

                if product_has_thioether:
                    # Check if this is a formation by seeing if any reactant lacks the pattern
                    formation_detected = False
                    for reactant in reactants:
                        if not checker.check_fg("Monosulfide", reactant):
                            formation_detected = True
                            break

                    if formation_detected:
                        has_thioether_formation = True
                        print(f"Detected thioether formation at depth {depth}")

            # Check for pyrimidinone formation
            if depth <= 2:
                # Check for pyrimidine ring with carbonyl group
                product_has_pyrimidine = checker.check_ring("pyrimidine", product)
                product_has_carbonyl = (
                    checker.check_fg("Ketone", product)
                    or checker.check_fg("Aldehyde", product)
                    or checker.check_fg("Carboxylic acid", product)
                    or checker.check_fg("Ester", product)
                    or checker.check_fg("Amide", product)
                )

                if product_has_pyrimidine and product_has_carbonyl:
                    # Check if this is a formation by seeing if any reactant lacks the pattern
                    formation_detected = False
                    for reactant in reactants:
                        reactant_has_pyrimidine = checker.check_ring("pyrimidine", reactant)
                        if not reactant_has_pyrimidine:
                            formation_detected = True
                            break

                    if formation_detected:
                        has_pyrimidinone_formation = True
                        print(f"Detected pyrimidinone formation at depth {depth}")

            # Check for convergent pattern at depth 1
            if depth == 1 and len(reactants) >= 2:
                has_convergent_pattern = True
                print(f"Detected convergent pattern at depth 1 with {len(reactants)} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the convergent synthesis with both key formations
    is_convergent_with_key_formations = (
        has_convergent_pattern
        and depth_2_fragments >= 2
        and (has_thioether_formation or has_pyrimidinone_formation)
    )

    print(f"Convergent: {has_convergent_pattern}, Fragments: {depth_2_fragments}")
    print(
        f"Thioether formation: {has_thioether_formation}, Pyrimidinone formation: {has_pyrimidinone_formation}"
    )
    print(
        f"Convergent synthesis with thioether and pyrimidinone: {is_convergent_with_key_formations}"
    )

    return is_convergent_with_key_formations
