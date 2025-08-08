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
    This function detects a strategy involving transformation of carbonyl groups to oxime ethers.
    """
    carbonyl_to_oxime_count = 0

    def dfs_traverse(node):
        nonlocal carbonyl_to_oxime_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for carbonyl to oxime transformation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if None not in reactant_mols and product_mol is not None:
                # Check if any reactant has C=O group
                carbonyl_pattern = Chem.MolFromSmarts("[C]=O")
                has_carbonyl = any(mol.HasSubstructMatch(carbonyl_pattern) for mol in reactant_mols)

                # Check if product has C=N-OCH3 group
                oxime_ether_pattern = Chem.MolFromSmarts("[C]=[N]-O[CH3]")
                has_oxime_ether = product_mol.HasSubstructMatch(oxime_ether_pattern)

                # Check if any reactant has methoxyamine
                methoxyamine_pattern = Chem.MolFromSmarts("[CH3]O[NH2]")
                has_methoxyamine = any(
                    mol.HasSubstructMatch(methoxyamine_pattern)
                    for mol in reactant_mols
                    if mol is not None
                )

                if has_carbonyl and has_oxime_ether and has_methoxyamine:
                    carbonyl_to_oxime_count += 1
                    print(f"Detected carbonyl to oxime transformation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if carbonyl to oxime transformations are detected
    result = carbonyl_to_oxime_count >= 1
    print(f"Total carbonyl to oxime transformations: {carbonyl_to_oxime_count}")
    return result
