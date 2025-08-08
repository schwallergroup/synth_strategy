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
    Detects if the synthesis route involves a decarbonylation step (removal of an aldehyde group).
    """
    decarbonylation_found = False

    def dfs_traverse(node):
        nonlocal decarbonylation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for decarbonylation
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and any(r for r in reactant_mols if r):
                # Check for aldehyde pattern in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=[O]")

                has_aldehyde = any(
                    r.HasSubstructMatch(aldehyde_pattern) for r in reactant_mols if r
                )

                # Check if product has fewer carbonyls than reactants
                if has_aldehyde:
                    # Count aldehyde groups in reactants and product
                    reactant_aldehyde_count = sum(
                        len(r.GetSubstructMatches(aldehyde_pattern)) for r in reactant_mols if r
                    )
                    product_aldehyde_count = len(product_mol.GetSubstructMatches(aldehyde_pattern))

                    if product_aldehyde_count < reactant_aldehyde_count:
                        decarbonylation_found = True
                        print("Decarbonylation found")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Decarbonylation strategy: {decarbonylation_found}")
    return decarbonylation_found
