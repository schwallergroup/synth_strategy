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
    This function detects preservation of bromine functional handle throughout the synthesis.
    It checks if the target molecule and a significant portion of intermediates contain a bromine atom.
    """
    # Track molecules with and without bromine
    molecules_with_bromine = 0
    total_non_starting_molecules = 0
    has_bromine_in_target = False

    def dfs_traverse(node, depth=0):
        nonlocal molecules_with_bromine, total_non_starting_molecules, has_bromine_in_target

        if node["type"] == "mol" and node["smiles"]:
            # Skip starting materials (in_stock=True)
            if not node.get("in_stock", False):
                mol_smiles = node["smiles"]
                total_non_starting_molecules += 1

                # Check if molecule contains bromine using RDKit
                mol = Chem.MolFromSmiles(mol_smiles)

                # Check for bromine atoms directly
                has_bromine = any(atom.GetSymbol() == "Br" for atom in mol.GetAtoms())

                # If molecule has bromine, increment counter
                if has_bromine:
                    molecules_with_bromine += 1
                    print(f"Molecule with bromine found at depth {depth}: {mol_smiles}")

                    # If this is the target molecule (depth 0), mark it
                    if depth == 0:
                        has_bromine_in_target = True
                else:
                    print(f"Molecule WITHOUT bromine found at depth {depth}: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Calculate threshold - at least 40% of intermediates should have bromine
    threshold_count = total_non_starting_molecules * 0.4

    # Check if target has bromine and enough intermediates have bromine
    bromine_preserved = has_bromine_in_target and (molecules_with_bromine >= threshold_count)

    print(
        f"Summary: {molecules_with_bromine}/{total_non_starting_molecules} non-starting molecules have bromine"
    )
    print(f"Target molecule has bromine: {has_bromine_in_target}")
    print(f"Threshold count: {threshold_count}")
    print(f"Bromine preservation strategy: {bromine_preserved}")

    return bromine_preserved
