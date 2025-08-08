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
    Detects a synthesis strategy involving diaryl ether formation.
    """
    has_diaryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_diaryl_ether_formation

        if node["type"] == "reaction":
            # Extract reactants and product from reaction
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(reactant_mols):
                # Check for diaryl ether formation
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")

                # Check if product has diaryl ether but reactants don't
                product_has_pattern = product_mol.HasSubstructMatch(diaryl_ether_pattern)
                reactants_have_pattern = any(
                    mol.HasSubstructMatch(diaryl_ether_pattern) for mol in reactant_mols
                )

                if product_has_pattern and not reactants_have_pattern:
                    has_diaryl_ether_formation = True
                    print(
                        f"Diaryl ether formation detected at depth {node.get('depth', 'unknown')}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if has_diaryl_ether_formation:
        print("Diaryl ether formation strategy detected")

    return has_diaryl_ether_formation
