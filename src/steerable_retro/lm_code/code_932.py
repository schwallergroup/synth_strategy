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
    This function detects a synthetic strategy involving multiple ether formations
    (at least 2) including benzyl and/or alkyl ethers.
    """
    ether_formation_count = 0

    def dfs_traverse(node):
        nonlocal ether_formation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ether formation
            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check if any reactant has a phenol or alcohol
                has_phenol_or_alcohol = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and (
                        reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][OH]"))
                        or reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C][OH]"))
                    ):
                        has_phenol_or_alcohol = True
                        break

                # Check if product has an ether
                has_ether = False
                if product_mol:
                    # Benzyl ether pattern
                    benzyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]")
                    # Alkyl ether pattern
                    alkyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][CH2]")

                    if product_mol.HasSubstructMatch(
                        benzyl_ether_pattern
                    ) or product_mol.HasSubstructMatch(alkyl_ether_pattern):
                        has_ether = True

                if has_phenol_or_alcohol and has_ether:
                    ether_formation_count += 1
                    print(f"Detected ether formation in reaction: {rsmi}")
            except:
                pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 ether formations are detected
    result = ether_formation_count >= 2
    print(f"Multiple ether formation strategy detected: {result} (count: {ether_formation_count})")
    return result
