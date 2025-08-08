#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects a strategy involving late-stage convergent coupling of two complex fragments,
    where at least one fragment contains a heterocycle.
    """
    found_late_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_convergent

        if (
            depth <= 1
            and node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Split reactants
            reactant_smiles = reactants_part.split(".")

            # Only consider reactions with at least 2 reactants
            if len(reactant_smiles) >= 2:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactant_smiles]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and all(r for r in reactant_mols):
                    # Check if at least one reactant has a heterocycle
                    heterocycle_patterns = [
                        Chem.MolFromSmarts("[n]1[c][n][c][c]1"),  # imidazole
                        Chem.MolFromSmarts("[s]1[c][n][c][c]1"),  # thiazole
                        Chem.MolFromSmarts("[o]1[c][c][c][c]1"),  # furan
                        Chem.MolFromSmarts("[nH]1[c][c][c][c]1"),  # pyrrole
                    ]

                    has_heterocycle = False
                    for r in reactant_mols:
                        if r and any(r.HasSubstructMatch(p) for p in heterocycle_patterns):
                            has_heterocycle = True
                            break

                    # Check if reactants are complex (at least 10 atoms)
                    complex_reactants = sum(1 for r in reactant_mols if r and r.GetNumAtoms() >= 10)

                    if has_heterocycle and complex_reactants >= 1:
                        found_late_convergent = True
                        print(f"Late-stage convergent coupling detected at depth {depth}")
                        print(f"Reactants: {reactant_smiles}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage convergent coupling strategy detected: {found_late_convergent}")

    return found_late_convergent
