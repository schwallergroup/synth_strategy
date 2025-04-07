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
    This function detects a synthetic strategy involving a pyridine scaffold with
    methoxy substituent that is maintained throughout the synthesis.
    """
    from rdkit import Chem

    pyridine_methoxy_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal pyridine_methoxy_count, total_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            total_reactions += 1
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for pyridine ring in the product
                has_pyridine = checker.check_ring("pyridine", product)

                if has_pyridine:
                    # Check specifically for methoxy group on pyridine
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        # Pattern for methoxy attached to aromatic carbon
                        methoxy_patt = Chem.MolFromSmarts("c-OC")

                        # First check if the molecule has a methoxy group
                        if mol.HasSubstructMatch(methoxy_patt):
                            # Then check if it's attached to the pyridine ring
                            pyridine_atoms = set()
                            try:
                                # Get all atoms in pyridine rings
                                matches = mol.GetSubstructMatches(Chem.MolFromSmarts("n1ccccc1"))
                                for match in matches:
                                    pyridine_atoms.update(match)

                                # Get all methoxy groups
                                methoxy_matches = mol.GetSubstructMatches(methoxy_patt)

                                # Check if any methoxy is attached to pyridine
                                for match in methoxy_matches:
                                    # match[0] is the aromatic carbon the methoxy is attached to
                                    if match[0] in pyridine_atoms:
                                        pyridine_methoxy_count += 1
                                        print(
                                            f"Pyridine with methoxy detected in product: {product}"
                                        )
                                        break
                            except Exception as e:
                                print(f"Error checking pyridine-methoxy attachment: {e}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if pyridine scaffold is maintained in most reactions (>80%)
    # and ensure there are enough reactions to make a meaningful assessment
    if total_reactions < 2:
        return False

    scaffold_maintained = pyridine_methoxy_count >= 0.8 * total_reactions

    print(f"Pyridine-methoxy count: {pyridine_methoxy_count}, Total reactions: {total_reactions}")
    print(f"Scaffold maintained: {scaffold_maintained}")

    return scaffold_maintained
