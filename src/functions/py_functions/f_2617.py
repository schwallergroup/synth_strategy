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


def main(route):
    """
    This function detects a synthetic strategy where a nitrogen heterocycle is formed
    in the late stage of the synthesis (low depth in the tree).
    """
    cyclization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclization_detected

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")
                    ]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        all(mol is not None for mol in reactants_mols)
                        and product_mol is not None
                    ):
                        # Count rings in reactants and product
                        reactants_ring_count = sum(
                            mol.GetRingInfo().NumRings() for mol in reactants_mols
                        )
                        product_ring_count = product_mol.GetRingInfo().NumRings()

                        # Check if a new ring is formed
                        if product_ring_count > reactants_ring_count:
                            # Check if the new ring contains nitrogen
                            nitrogen_pattern = Chem.MolFromSmarts("[#7]")
                            if product_mol.HasSubstructMatch(nitrogen_pattern):
                                # Check if the nitrogen is part of a ring
                                for atom in product_mol.GetAtoms():
                                    if atom.GetAtomicNum() == 7 and atom.IsInRing():
                                        cyclization_detected = True
                                        print(
                                            f"Late-stage nitrogen heterocycle formation detected at depth {depth}"
                                        )
                                        break
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return cyclization_detected
