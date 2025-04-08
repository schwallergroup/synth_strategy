#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects benzisoxazole formation via oxime intermediate in the synthesis route.
    Looks for the sequence: ketone -> oxime -> benzisoxazole
    """
    # Track molecules and reactions in the synthesis path
    synthesis_paths = []
    current_path = []

    def dfs_traverse(node, depth=0):
        # Add current node to path
        if node["type"] == "mol":
            mol_info = {"type": "mol", "smiles": node["smiles"], "depth": depth}
            current_path.append(mol_info)
        elif node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            rxn_info = {
                "type": "reaction",
                "rsmi": rsmi,
                "reactants": reactants,
                "product": product,
                "depth": depth,
            }
            current_path.append(rxn_info)

        # If this is a leaf node (starting material), save the current path
        if node["type"] == "mol" and node.get("in_stock", False):
            synthesis_paths.append(list(current_path))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

        # Remove current node from path when backtracking
        if len(current_path) > 0:
            current_path.pop()

    # Start traversal
    dfs_traverse(route)

    # Check each synthesis path for the sequence
    for path in synthesis_paths:
        # Sort by depth to get chronological synthesis order (high to low depth)
        path.sort(key=lambda x: x["depth"], reverse=True)

        # Look for the sequence: ketone -> oxime -> benzisoxazole
        ketone_steps = []
        oxime_steps = []
        benzisoxazole_steps = []

        for i, step in enumerate(path):
            if step["type"] == "mol":
                mol_smiles = step["smiles"]

                # Check for ketone
                if checker.check_fg("Ketone", mol_smiles):
                    print(f"Found ketone at step {i}: {mol_smiles}")
                    ketone_steps.append((i, mol_smiles))

                # Check for oxime
                if checker.check_fg("Oxime", mol_smiles):
                    print(f"Found oxime at step {i}: {mol_smiles}")
                    oxime_steps.append((i, mol_smiles))

                # Check for benzisoxazole - isoxazole fused to benzene
                if checker.check_ring("isoxazole", mol_smiles) and "c1ccc" in mol_smiles:
                    print(f"Found benzisoxazole at step {i}: {mol_smiles}")
                    benzisoxazole_steps.append((i, mol_smiles))

            elif step["type"] == "reaction":
                # Check for oxime formation reaction
                if any(
                    checker.check_fg("Ketone", r) for r in step["reactants"]
                ) and checker.check_fg("Oxime", step["product"]):
                    print(f"Found oxime formation reaction: {step['rsmi']}")

                # Check for benzisoxazole formation reaction
                if (
                    any(checker.check_fg("Oxime", r) for r in step["reactants"])
                    and checker.check_ring("isoxazole", step["product"])
                    and "c1ccc" in step["product"]
                ):
                    print(f"Found benzisoxazole formation reaction: {step['rsmi']}")

        # Check if we have all three components
        if ketone_steps and oxime_steps and benzisoxazole_steps:
            # Check for valid sequences (ketone -> oxime -> benzisoxazole)
            for k_idx, k_smiles in ketone_steps:
                for o_idx, o_smiles in oxime_steps:
                    for b_idx, b_smiles in benzisoxazole_steps:
                        # Check correct sequence order
                        if k_idx < o_idx < b_idx:
                            print(f"Found sequence in correct order: {k_idx} -> {o_idx} -> {b_idx}")

                            # Check for structural similarity between steps
                            # This is a simplified check - in a real implementation,
                            # you would use atom mapping to ensure the same fragment is transformed
                            k_mol = Chem.MolFromSmiles(k_smiles)
                            o_mol = Chem.MolFromSmiles(o_smiles)
                            b_mol = Chem.MolFromSmiles(b_smiles)

                            if k_mol and o_mol and b_mol:
                                # Check if the oxime has similar structure to the ketone
                                # and the benzisoxazole has similar structure to the oxime
                                # This is a very basic check and could be improved
                                k_frags = k_smiles.split("C(=O)")
                                o_frags = o_smiles.split("C(=NO)")

                                if len(k_frags) > 1 and len(o_frags) > 1:
                                    # If fragments around the ketone/oxime are similar
                                    if any(kf in o_smiles for kf in k_frags if kf):
                                        print("Found complete sequence with structural continuity")
                                        return True

    return False
