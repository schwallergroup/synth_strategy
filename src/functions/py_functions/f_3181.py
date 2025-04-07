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
    Detects if the synthesis route involves heterocycle formation.
    """
    has_heterocycle_formation = False

    def dfs_traverse(node):
        nonlocal has_heterocycle_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycle formation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol is not None:
                # Count rings in reactants and product
                reactant_ring_info = [
                    Chem.GetSSSR(mol) for mol in reactant_mols if mol is not None
                ]
                reactant_ring_count = sum(len(rings) for rings in reactant_ring_info)

                product_ring_info = Chem.GetSSSR(product_mol)
                product_ring_count = len(product_ring_info)

                # Check if product has more rings than reactants combined
                if product_ring_count > reactant_ring_count:
                    # Check if new ring contains heteroatom
                    for ring in product_ring_info:
                        ring_atoms = list(ring)
                        ring_smiles = Chem.MolFragmentToSmiles(product_mol, ring_atoms)
                        if any(
                            atom in ring_smiles
                            for atom in ["n", "N", "o", "O", "s", "S"]
                        ):
                            has_heterocycle_formation = True
                            print(f"Detected heterocycle formation: {ring_smiles}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_heterocycle_formation
