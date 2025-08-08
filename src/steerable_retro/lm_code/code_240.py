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
    Detects if the synthesis involves a Wittig-type olefination connecting two aromatic rings
    """
    wittig_detected = False

    def dfs_traverse(node):
        nonlocal wittig_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phosphonate pattern in reactants
            phosphonate_pattern = Chem.MolFromSmarts("[P](=[O])([O])[CH2][c]")
            # Check for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[CH]=[O]")

            # Check for olefin connecting two aromatic rings in product
            olefin_pattern = Chem.MolFromSmarts("[c]/[CH]=[CH]/[c]")

            has_phosphonate = False
            has_aldehyde = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phosphonate_pattern):
                        has_phosphonate = True
                    if mol and mol.HasSubstructMatch(aldehyde_pattern):
                        has_aldehyde = True
                except:
                    continue

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if (
                    has_phosphonate
                    and has_aldehyde
                    and prod_mol
                    and prod_mol.HasSubstructMatch(olefin_pattern)
                ):
                    print("Wittig-type olefination detected")
                    wittig_detected = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return wittig_detected
