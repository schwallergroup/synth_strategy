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
    This function detects enaminone formation from a ketone and DMF acetal.
    """
    enaminone_formation_detected = False

    def dfs_traverse(node):
        nonlocal enaminone_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant is a ketone
            ketone_pattern = Chem.MolFromSmarts("[C][C](=O)[c,C]")

            # Check if another reactant is DMF acetal or similar
            dmf_acetal_pattern = Chem.MolFromSmarts("CO[CH](OC)N(C)C")

            # Check if product is an enaminone
            enaminone_pattern = Chem.MolFromSmarts("CN(C)C=CC(=O)[c,C]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(mol and mol.HasSubstructMatch(ketone_pattern) for mol in reactant_mols)
                and any(mol and mol.HasSubstructMatch(dmf_acetal_pattern) for mol in reactant_mols)
                and product_mol.HasSubstructMatch(enaminone_pattern)
            ):
                enaminone_formation_detected = True
                print("Detected enaminone formation from ketone and DMF acetal")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return enaminone_formation_detected
