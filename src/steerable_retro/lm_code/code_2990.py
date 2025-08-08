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
    This function detects if the synthesis involves sulfonamide formation from sulfonyl chloride.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains sulfonamide
                product_mol = Chem.MolFromSmiles(product)
                sulfonamide_pattern = Chem.MolFromSmarts("[c][S](=[O])(=[O])[N]")

                # Check if reactants contain sulfonyl chloride and amine
                has_sulfonyl_chloride = False
                has_amine = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        sulfonyl_chloride_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[Cl]")
                        amine_pattern = Chem.MolFromSmarts("[N]")

                        if reactant_mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                            has_sulfonyl_chloride = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                # If product has sulfonamide and reactants have sulfonyl chloride and amine
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(sulfonamide_pattern)
                    and has_sulfonyl_chloride
                    and has_amine
                ):
                    has_sulfonamide_formation = True
                    print(f"Detected sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_sulfonamide_formation
