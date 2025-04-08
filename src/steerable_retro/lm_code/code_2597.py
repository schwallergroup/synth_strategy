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
    This function detects if a heterocyclic ring formation occurs in the late stage
    (final 2 steps) of the synthesis.
    """
    late_stage_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_ring_formation

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if all(mol is not None for mol in reactant_mols) and product_mol is not None:
                # Count rings in reactants and product
                reactant_ring_count = sum(
                    rdMolDescriptors.CalcNumRings(mol) for mol in reactant_mols
                )
                product_ring_count = rdMolDescriptors.CalcNumRings(product_mol)

                # Check if product has more rings than reactants
                if product_ring_count > reactant_ring_count:
                    # Check if any new ring contains a heteroatom
                    n_pattern = Chem.MolFromSmarts("[#7]")  # Nitrogen
                    o_pattern = Chem.MolFromSmarts("[#8]")  # Oxygen
                    s_pattern = Chem.MolFromSmarts("[#16]")  # Sulfur

                    has_n = product_mol.HasSubstructMatch(n_pattern)
                    has_o = product_mol.HasSubstructMatch(o_pattern)
                    has_s = product_mol.HasSubstructMatch(s_pattern)

                    if has_n or has_o or has_s:
                        print(f"Detected heterocyclic ring formation at depth {depth}")
                        late_stage_ring_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_ring_formation
