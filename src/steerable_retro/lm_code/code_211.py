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
    This function detects if the synthesis route includes a sulfonamide formation step
    from an amine and sulfonyl chloride.
    """
    sulfonamide_formation_found = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                product_mol = Chem.MolFromSmiles(product)

                # Check for sulfonamide pattern in product
                sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=[O])(=[O])")
                sulfonyl_chloride_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[Cl]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                if (
                    product_mol
                    and sulfonamide_pattern
                    and product_mol.HasSubstructMatch(sulfonamide_pattern)
                ):
                    # Check if reactants include sulfonyl chloride and amine
                    has_sulfonyl_chloride = False
                    has_amine = False

                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            if sulfonyl_chloride_pattern and r_mol.HasSubstructMatch(
                                sulfonyl_chloride_pattern
                            ):
                                has_sulfonyl_chloride = True
                            if amine_pattern and r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_sulfonyl_chloride and has_amine:
                        sulfonamide_formation_found = True
                        print("Sulfonamide formation detected")

            except Exception as e:
                print(f"Error in sulfonamide detection: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_formation_found
