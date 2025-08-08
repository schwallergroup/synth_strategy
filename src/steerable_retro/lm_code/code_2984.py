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
    Detects if the synthesis route involves an olefination reaction
    (formation of C=C bond using an aldehyde).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth > 1:  # Early stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde in reactants
                aldehyde_found = False
                for reactant in reactants:
                    try:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                        if reactant_mol and reactant_mol.HasSubstructMatch(aldehyde_pattern):
                            aldehyde_found = True
                            break
                    except:
                        continue

                # Check for alkene in product
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    alkene_pattern = Chem.MolFromSmarts("C=C")
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(alkene_pattern)
                        and aldehyde_found
                    ):
                        print(f"Detected olefination reaction at depth {depth}")
                        result = True
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
