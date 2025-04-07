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
    This function detects if the synthesis route involves formation of a heterocyclic ring system
    (particularly pyridine) in the early stages of the synthesis.
    """
    heterocycle_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and depth >= 4:  # Early stage (high depth)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count rings in reactants and product
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and all(reactant_mols):
                    # Count nitrogen-containing rings in reactants
                    reactant_n_rings = 0
                    for mol in reactant_mols:
                        for ring in Chem.GetSSSR(mol):
                            ring_atoms = set(ring)
                            for atom_idx in ring_atoms:
                                atom = mol.GetAtomWithIdx(atom_idx)
                                if atom.GetAtomicNum() == 7:  # Nitrogen
                                    reactant_n_rings += 1
                                    break

                    # Count nitrogen-containing rings in product
                    product_n_rings = 0
                    for ring in Chem.GetSSSR(product_mol):
                        ring_atoms = set(ring)
                        for atom_idx in ring_atoms:
                            atom = product_mol.GetAtomWithIdx(atom_idx)
                            if atom.GetAtomicNum() == 7:  # Nitrogen
                                product_n_rings += 1
                                break

                    # If product has more nitrogen rings than reactants combined, it's a heterocycle formation
                    if product_n_rings > reactant_n_rings:
                        print(f"Heterocycle formation detected at depth {depth}")
                        heterocycle_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_formation_detected
