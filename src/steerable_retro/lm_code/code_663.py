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
    This function detects if the synthesis preserves stereochemistry throughout
    (no stereocenter modifications).
    """
    # This is a complex detection that would require tracking stereocenters through reactions
    # For simplicity, we'll implement a basic version that checks if stereocenters are present
    # in both reactants and products for each reaction

    preserves_stereo = True

    def dfs_traverse(node):
        nonlocal preserves_stereo

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                # Check stereocenters in reactants
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                reactants_stereo = False
                if reactants_mol:
                    for atom in reactants_mol.GetAtoms():
                        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                            reactants_stereo = True
                            break

                # Check stereocenters in products
                products_mol = Chem.MolFromSmiles(products_smiles)
                products_stereo = False
                if products_mol:
                    for atom in products_mol.GetAtoms():
                        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                            products_stereo = True
                            break

                # If stereocenters are present in one but not the other, stereochemistry is not preserved
                if reactants_stereo != products_stereo:
                    preserves_stereo = False
                    print("Found stereochemistry modification")
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Preserves stereochemistry: {preserves_stereo}")
    return preserves_stereo
