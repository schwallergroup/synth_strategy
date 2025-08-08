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
    This function detects if the synthesis route starts with a hydrazine condensation.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth >= 2:  # Early stage (high depth)
            # Check if this is a hydrazine condensation
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                try:
                    # Check for hydrazine pattern in reactants
                    for reactant in reactants_part.split("."):
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            hydrazine_pattern = Chem.MolFromSmarts("[N]-[N]")
                            if mol.HasSubstructMatch(hydrazine_pattern):
                                # Check if product has C=N-N pattern (hydrazone)
                                product_mol = Chem.MolFromSmiles(product_part)
                                if product_mol:
                                    hydrazone_pattern = Chem.MolFromSmarts("[C]=[N]-[N]")
                                    if product_mol.HasSubstructMatch(hydrazone_pattern):
                                        result = True
                                        print(
                                            f"Found early hydrazine condensation at depth {depth}"
                                        )
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return result
