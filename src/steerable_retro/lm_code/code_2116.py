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
    This function detects if the synthetic route involves protection and/or deprotection of a phenol group.
    Specifically looking for MOM (methoxymethyl) protection.
    """
    phenol_pattern = Chem.MolFromSmarts("Oc1ccccc1")
    protected_phenol_pattern = Chem.MolFromSmarts("COCOc1ccccc1")
    has_phenol_protection = False

    def dfs_traverse(node):
        nonlocal has_phenol_protection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check for protection (phenol -> protected phenol)
                    if (
                        reactants_mol
                        and product_mol
                        and reactants_mol.HasSubstructMatch(phenol_pattern)
                        and product_mol.HasSubstructMatch(protected_phenol_pattern)
                    ):
                        print(f"Detected phenol protection in reaction: {rsmi}")
                        has_phenol_protection = True

                    # Check for deprotection (protected phenol -> phenol)
                    if (
                        reactants_mol
                        and product_mol
                        and reactants_mol.HasSubstructMatch(protected_phenol_pattern)
                        and product_mol.HasSubstructMatch(phenol_pattern)
                    ):
                        print(f"Detected phenol deprotection in reaction: {rsmi}")
                        has_phenol_protection = True
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_phenol_protection
