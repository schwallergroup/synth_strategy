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
    This function detects pyrazole ring formation via condensation of hydrazine with an enaminone.
    """
    pyrazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal pyrazole_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if hydrazine is a reactant
            hydrazine_pattern = Chem.MolFromSmarts("[NH2][NH2]")

            # Check if product contains pyrazole
            pyrazole_pattern = Chem.MolFromSmarts("c1[nH]ncc1")

            # Check if one reactant contains enaminone pattern
            enaminone_pattern = Chem.MolFromSmarts("CN(C)C=CC(=O)c")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(mol and mol.HasSubstructMatch(hydrazine_pattern) for mol in reactant_mols)
                and any(mol and mol.HasSubstructMatch(enaminone_pattern) for mol in reactant_mols)
                and product_mol.HasSubstructMatch(pyrazole_pattern)
            ):
                pyrazole_formation_detected = True
                print("Detected pyrazole formation from enaminone and hydrazine")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrazole_formation_detected
