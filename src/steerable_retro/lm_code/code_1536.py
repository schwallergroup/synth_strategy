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
    This function detects SNAr reaction forming aryl ether using fluoronitrobenzene.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for fluoronitrobenzene in reactants
            fluoro_nitro_pattern = Chem.MolFromSmarts("c1(F)ccc(cc1)[N+](=O)[O-]")
            phenol_pattern = Chem.MolFromSmarts("c1(O)ccccc1")
            aryl_ether_pattern = Chem.MolFromSmarts("c1(Oc2ccccc2)ccc(cc1)[N+](=O)[O-]")

            fluoro_nitro_present = False
            phenol_present = False

            for reactant in reactants_smiles:
                r_mol = Chem.MolFromSmiles(reactant)
                if r_mol:
                    if r_mol.HasSubstructMatch(fluoro_nitro_pattern):
                        fluoro_nitro_present = True
                    if r_mol.HasSubstructMatch(phenol_pattern):
                        phenol_present = True

            # Check for aryl ether in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern):
                if fluoro_nitro_present and phenol_present:
                    snar_detected = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"SNAr aryl ether formation detected: {snar_detected}")
    return snar_detected
