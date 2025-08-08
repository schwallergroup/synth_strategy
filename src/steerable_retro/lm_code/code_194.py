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
    This function detects a strategy involving a β-lactam core with phenyl substituent
    that is maintained throughout the synthesis.
    """
    # Track if target molecule and at least one intermediate have β-lactam with phenyl
    found_beta_lactam_with_phenyl = False

    def is_beta_lactam_with_phenyl(smiles):
        """Check if a molecule contains a β-lactam with phenyl substituent"""
        try:
            # β-lactam is a 4-membered ring with nitrogen and carbonyl
            # We need to check for both the ring structure and the phenyl group
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # Check for β-lactam structure (4-membered ring with N and C=O)
            beta_lactam_pattern = Chem.MolFromSmarts("[#6]1[#6](=[#8])[#7][#6]1")
            has_beta_lactam = mol.HasSubstructMatch(beta_lactam_pattern)

            # Check for phenyl group
            phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
            has_phenyl = mol.HasSubstructMatch(phenyl_pattern)

            # Check if the phenyl is connected to the β-lactam
            if has_beta_lactam and has_phenyl:
                # Get the atoms in the β-lactam
                beta_lactam_matches = mol.GetSubstructMatches(beta_lactam_pattern)
                if not beta_lactam_matches:
                    return False

                # Get the atoms in the phenyl group
                phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)
                if not phenyl_matches:
                    return False

                # Check if any atom in the β-lactam is bonded to any atom in the phenyl
                for beta_match in beta_lactam_matches:
                    for phenyl_match in phenyl_matches:
                        for beta_atom_idx in beta_match:
                            beta_atom = mol.GetAtomWithIdx(beta_atom_idx)
                            for neighbor in beta_atom.GetNeighbors():
                                if neighbor.GetIdx() in phenyl_match:
                                    print(f"Found β-lactam with phenyl in: {smiles}")
                                    return True

            return False
        except Exception as e:
            print(f"Error checking for β-lactam with phenyl in {smiles}: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal found_beta_lactam_with_phenyl

        if node["type"] == "mol" and "smiles" in node:
            # Skip checking for starting materials
            if not node.get("in_stock", False):
                if is_beta_lactam_with_phenyl(node["smiles"]):
                    # If this is the target molecule (depth 0), mark as found
                    if depth == 0:
                        found_beta_lactam_with_phenyl = True
                else:
                    # If any non-starting material doesn't have the pattern, print it
                    print(f"Molecule without β-lactam-phenyl found: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Target molecule has β-lactam with phenyl: {found_beta_lactam_with_phenyl}")

    return found_beta_lactam_with_phenyl
