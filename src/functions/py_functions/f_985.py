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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if a bromopyridine moiety is preserved throughout the synthesis.
    It checks if the target molecule contains a bromopyridine and if this structure
    is preserved in all non-starting material molecules and through all reactions.
    """
    # Helper function to check if a molecule contains a bromopyridine moiety
    def has_bromopyridine(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                print(f"Could not parse SMILES: {smiles}")
                return False

            # Check if molecule has pyridine ring
            has_pyridine = checker.check_ring("pyridine", smiles)
            if not has_pyridine:
                return False

            # Check if molecule has bromine attached to pyridine
            # Get pyridine atom indices
            pyridine_indices = checker.get_ring_atom_indices("pyridine", smiles)
            if not pyridine_indices:
                return False

            # Check if any bromine is attached to pyridine atoms
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 35:  # Bromine
                    # Check if bromine is attached to any pyridine atom
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        for ring_indices in pyridine_indices:
                            if neighbor_idx in ring_indices:
                                print(f"Found bromopyridine in: {smiles}")
                                return True
            return False
        except Exception as e:
            print(f"Error checking bromopyridine: {e}")
            return False

    # Check if the target molecule has a bromopyridine
    target_has_bromopyridine = False
    if route["type"] == "mol" and "smiles" in route:
        target_has_bromopyridine = has_bromopyridine(route["smiles"])

    if not target_has_bromopyridine:
        print(
            f"Target molecule does not contain bromopyridine: {route.get('smiles', '')}"
        )
        return False

    # Track if bromopyridine is preserved throughout synthesis
    preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal preserved

        # Skip further checks if we already know it's not preserved
        if not preserved:
            return

        if node["type"] == "mol":
            # Skip starting materials
            if node.get("in_stock", False):
                print(f"Skipping starting material: {node.get('smiles', '')}")
                return

            # Check if non-starting material molecules have bromopyridine
            if "smiles" in node and not has_bromopyridine(node["smiles"]):
                print(f"Molecule without bromopyridine found: {node['smiles']}")
                preserved = False
                return

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # For reaction nodes, check if bromopyridine is preserved
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if any reactant has bromopyridine
                reactant_has_bromopyridine = any(
                    has_bromopyridine(r) for r in reactants
                )

                # If any reactant has bromopyridine, product should have it too
                if reactant_has_bromopyridine and not has_bromopyridine(product):
                    print(f"Bromopyridine lost in reaction: {rsmi}")
                    preserved = False
                    return

                # If product has bromopyridine but no reactant does, it was created
                if has_bromopyridine(product) and not reactant_has_bromopyridine:
                    print(f"Bromopyridine created in reaction: {rsmi}")
                    preserved = False
                    return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print(f"Starting traversal of synthesis route")
    dfs_traverse(route)

    print(f"Bromopyridine preserved throughout synthesis: {preserved}")
    return preserved
