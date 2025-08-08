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
    Detects if the synthesis involves formation of a sulfonamide bond.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Sulfonamide pattern
            sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product = Chem.MolFromSmiles(product_part)

                if product and sulfonamide_pattern:
                    # Check if sulfonamide is in product but not in all reactants
                    product_has_pattern = product.HasSubstructMatch(sulfonamide_pattern)
                    reactants_have_pattern = all(
                        r and r.HasSubstructMatch(sulfonamide_pattern) for r in reactants if r
                    )

                    if product_has_pattern and not reactants_have_pattern:
                        print("Detected sulfonamide formation")
                        has_sulfonamide_formation = True
            except:
                print("Error in SMILES processing for sulfonamide detection")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_sulfonamide_formation
