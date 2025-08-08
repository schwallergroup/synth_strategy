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
    This function detects if the synthetic route employs thiazole formation
    from α-bromoketone and thiourea.
    """
    thiazole_formed = False

    def dfs_traverse(node):
        nonlocal thiazole_formed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Define SMARTS patterns for α-bromoketone, thiourea, and thiazole
            bromoketone_pattern = Chem.MolFromSmarts("C(=O)C(Br)")
            thiourea_pattern = Chem.MolFromSmarts("NC(=S)N")
            thiazole_pattern = Chem.MolFromSmarts("c1nc(N)sc1")

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check if reactants contain α-bromoketone and thiourea, and product contains thiazole
                    if (
                        reactants_mol.HasSubstructMatch(bromoketone_pattern)
                        and reactants_mol.HasSubstructMatch(thiourea_pattern)
                        and product_mol.HasSubstructMatch(thiazole_pattern)
                    ):
                        print("Detected thiazole formation from α-bromoketone and thiourea")
                        thiazole_formed = True
            except:
                print("Error processing SMILES in thiazole formation detection")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return thiazole_formed
