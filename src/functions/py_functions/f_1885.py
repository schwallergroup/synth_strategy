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
    Detects if aromatic chlorides are maintained throughout the synthesis
    """
    # Track molecules with aromatic chlorides at each depth
    molecules_with_chlorides = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule contains aromatic chlorides
            mol_smiles = node["smiles"]
            if checker.check_fg("Aromatic halide", mol_smiles) and "Cl" in mol_smiles:
                # Count aromatic chlorides specifically
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    ar_cl_pattern = Chem.MolFromSmarts("c-[Cl]")
                    count = len(mol.GetSubstructMatches(ar_cl_pattern))
                    if count > 0:
                        print(
                            f"Found {count} aromatic chlorides at depth {depth} in molecule: {mol_smiles}"
                        )
                        molecules_with_chlorides.append((depth, count, mol_smiles))

        # Continue traversing the synthesis tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root (target molecule)
    dfs_traverse(route)

    # If we found at least one molecule with aromatic chlorides
    if molecules_with_chlorides:
        # Sort by depth (lowest depth = target molecule, highest depth = starting materials)
        molecules_with_chlorides.sort(key=lambda x: x[0])

        # Get the target molecule (lowest depth)
        target_depth, target_count, target_smiles = molecules_with_chlorides[0]

        # Check if target molecule has aromatic chlorides
        if target_count > 0:
            print(f"Target molecule has {target_count} aromatic chlorides")

            # Check if there are starting materials or intermediates with aromatic chlorides
            if len(molecules_with_chlorides) > 1:
                # Check if any starting material has more aromatic chlorides than the target
                for depth, count, smiles in molecules_with_chlorides[1:]:
                    if count > target_count:
                        print(
                            f"Molecule at depth {depth} has {count} aromatic chlorides, which exceeds target's {target_count}"
                        )
                        return False

                # If we get here, all molecules have equal or fewer aromatic chlorides than the target
                print(
                    f"Aromatic chlorides are maintained: target has {target_count}, no intermediate or starting material exceeds this"
                )
                return True

    # Default case: no aromatic chlorides maintained
    return False
