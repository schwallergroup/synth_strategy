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
    This function detects a strategy involving multiple ether formation reactions
    (C-O-C bonds) at different stages of the synthesis.
    """
    ether_formation_count = 0
    ether_formation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_count, ether_formation_depths

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ether formation (C-O-C bond that wasn't present in reactants)
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactants_mols):
                # Look for ether pattern in product
                ether_pattern = Chem.MolFromSmarts("[#6][#8][#6]")
                product_ethers = product_mol.GetSubstructMatches(ether_pattern)

                # Check each reactant for the same ether bonds
                new_ether_formed = False
                for ether_match in product_ethers:
                    ether_atoms = set(ether_match)

                    # Check if this ether bond exists in any reactant
                    exists_in_reactants = False
                    for r_mol in reactants_mols:
                        if r_mol.HasSubstructMatch(ether_pattern):
                            for r_match in r_mol.GetSubstructMatches(ether_pattern):
                                if set(r_match) == ether_atoms:
                                    exists_in_reactants = True
                                    break

                    if not exists_in_reactants:
                        new_ether_formed = True
                        break

                if new_ether_formed:
                    ether_formation_count += 1
                    ether_formation_depths.append(depth)
                    print(f"Ether formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy criteria: At least 2 ether formations at different depths
    result = ether_formation_count >= 2 and len(set(ether_formation_depths)) >= 2
    print(f"Sequential ether formations strategy detected: {result}")
    print(f"Total ether formations: {ether_formation_count} at depths {ether_formation_depths}")

    return result
