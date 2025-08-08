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
    Detects a late-stage fragment coupling via N-alkylation between two aromatic systems.
    """
    found_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions (low depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if len(reactants) >= 2:  # Need at least two reactants for coupling
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if product_mol:
                        # Check if reactants are aromatic
                        aromatic_reactants = 0
                        for mol in reactant_mols:
                            if mol and mol.GetNumAtoms() > 6:  # Non-trivial fragment
                                aromatic_atoms = sum(
                                    1 for atom in mol.GetAtoms() if atom.GetIsAromatic()
                                )
                                if aromatic_atoms >= 6:  # At least a benzene ring
                                    aromatic_reactants += 1

                        # Check for N-alkylation pattern
                        if aromatic_reactants >= 2:
                            # Look for patterns indicating N-alkylation
                            n_alkylation_pattern = Chem.MolFromSmarts("[c][N][C][c]")
                            if product_mol.HasSubstructMatch(n_alkylation_pattern):
                                # Check if this bond is newly formed
                                found_n_alkylation = True
                                print(f"Found late-stage aromatic N-alkylation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_n_alkylation
