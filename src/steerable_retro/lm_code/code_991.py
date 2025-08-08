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
    This function detects aromatic halogenation reactions.
    """
    aromatic_halogenation_detected = False

    def dfs_traverse(node):
        nonlocal aromatic_halogenation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for halogen sources in reactants
            halogen_source_patterns = [
                Chem.MolFromSmarts("II"),  # I2
                Chem.MolFromSmarts("BrBr"),  # Br2
                Chem.MolFromSmarts("ClCl"),  # Cl2
            ]

            has_halogen_source = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        for pattern in halogen_source_patterns:
                            if pattern and mol.HasSubstructMatch(pattern):
                                has_halogen_source = True
                                break
                except:
                    continue

            # Check for aromatic-halogen bond in product
            aromatic_halogen_pattern = Chem.MolFromSmarts("c-[Br,I,Cl]")

            has_aromatic_halogen = False
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(aromatic_halogen_pattern):
                    has_aromatic_halogen = True
            except:
                pass

            if has_halogen_source and has_aromatic_halogen:
                print("Detected aromatic halogenation")
                aromatic_halogenation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return aromatic_halogenation_detected
