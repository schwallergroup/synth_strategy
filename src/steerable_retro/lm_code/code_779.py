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
    This function detects if an early step in the synthesis involves
    halogenation of an aromatic ring.
    """
    # Track depth of aromatic halogenation
    halogenation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]

                    if product_mol and reactant_mols:
                        # Check for aromatic halogenation
                        aromatic_halogen_pattern = Chem.MolFromSmarts("[c][Br,Cl,I,F]")

                        # Check if product has aromatic halogen that wasn't in reactants
                        if product_mol.HasSubstructMatch(aromatic_halogen_pattern):
                            # Check if this is a new halogenation
                            reactant_has_pattern = any(
                                mol.HasSubstructMatch(aromatic_halogen_pattern)
                                for mol in reactant_mols
                            )

                            if not reactant_has_pattern:
                                halogenation_depth = depth
                except:
                    print("Error processing molecule in early_aromatic_halogenation")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if halogenation is in the first 30% of steps
    if halogenation_depth is not None and halogenation_depth >= max_depth * 0.7:
        print(f"Detected early aromatic halogenation at depth {halogenation_depth} of {max_depth}")
        return True
    return False
