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
    This function detects if the synthesis involves multiple chloride installations
    as leaving groups (e.g., alcohol to chloride, nitro to chloride).
    """
    chloride_installations = 0

    def dfs_traverse(node):
        nonlocal chloride_installations

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for alcohol to chloride transformation
            alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")
            nitro_pattern = Chem.MolFromSmarts("[#6][N+](=[O])[O-]")

            alcohol_present = False
            nitro_present = False

            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol:
                        if r_mol.HasSubstructMatch(alcohol_pattern):
                            alcohol_present = True
                        if r_mol.HasSubstructMatch(nitro_pattern):
                            nitro_present = True
                except:
                    continue

            # Check if product has chloride where alcohol or nitro was
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                chloride_pattern = Chem.MolFromSmarts("[#6]Cl")

                if p_mol and p_mol.HasSubstructMatch(chloride_pattern):
                    if alcohol_present or nitro_present:
                        chloride_installations += 1
                        print(f"Detected chloride installation (count: {chloride_installations})")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return chloride_installations >= 2  # Return True if at least 2 chloride installations
