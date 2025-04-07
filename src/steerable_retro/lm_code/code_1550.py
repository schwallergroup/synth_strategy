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
    Detects synthesis routes that build a piperazine scaffold with sequential N-alkylations.
    """
    # Track if we found piperazine and alkylation reactions
    found_piperazine = False
    alkylation_count = 0

    def dfs_traverse(node):
        nonlocal found_piperazine, alkylation_count

        if node["type"] == "mol":
            # Check if molecule contains piperazine
            if node.get("smiles"):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
                    if mol.HasSubstructMatch(piperazine_pattern):
                        found_piperazine = True
                        print("Found piperazine scaffold")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alkylation reaction (N-alkylation)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Look for alkyl halide pattern in reactants
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6]-[Br,Cl,I]")
                amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O)]")

                has_alkyl_halide = any(
                    mol.HasSubstructMatch(alkyl_halide_pattern) for mol in reactant_mols
                )
                has_amine = any(mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)

                if has_alkyl_halide and has_amine:
                    alkylation_count += 1
                    print(f"Found N-alkylation reaction, total count: {alkylation_count}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found piperazine and at least 2 alkylation reactions
    result = found_piperazine and alkylation_count >= 2
    print(f"Piperazine scaffold with sequential alkylations strategy detected: {result}")
    return result
