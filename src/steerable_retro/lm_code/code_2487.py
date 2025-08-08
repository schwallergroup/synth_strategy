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
    This function detects if the synthetic route involves a convergent synthesis
    where two complex fragments (>10 atoms each) are joined.
    """
    # Flag to track if we found the pattern
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        # Only process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and products
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Parse reactants and product
            reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
            product = Chem.MolFromSmiles(product_str)

            if product and len(reactants) >= 2 and all(r for r in reactants):
                # Count atoms in each reactant
                reactant_atom_counts = [r.GetNumAtoms() for r in reactants if r]

                # Check if we have at least two complex fragments (>10 atoms each)
                complex_fragments = [count for count in reactant_atom_counts if count > 10]

                if len(complex_fragments) >= 2:
                    print(f"Found convergent synthesis with fragment sizes: {reactant_atom_counts}")
                    found_pattern = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
