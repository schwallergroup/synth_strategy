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
    This function detects if the synthesis involves formation of a heterocyclic system,
    particularly focusing on reactions that increase the ring count.
    """
    heterocycle_formation_detected = False

    def dfs_traverse(node):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Count rings in reactants
                    reactant_ring_count = 0
                    for r_smiles in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol:
                            reactant_ring_count += r_mol.GetRingInfo().NumRings()

                    # Count rings in product
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol:
                        product_ring_count = p_mol.GetRingInfo().NumRings()

                        # Check if product has more rings than reactants combined
                        if product_ring_count > reactant_ring_count:
                            # Check if any of the new rings contains a heteroatom
                            for atom in p_mol.GetAtoms():
                                if atom.IsInRing() and atom.GetAtomicNum() not in [
                                    1,
                                    6,
                                ]:  # Not H or C
                                    heterocycle_formation_detected = True
                                    print(f"Heterocycle formation detected: {rsmi}")
                                    break
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_formation_detected
