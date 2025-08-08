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
    This function detects a silyl protection-deprotection sequence in the synthesis.
    Specifically looking for TBS (tert-butyldimethylsilyl) protection of alcohols.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                products = [Chem.MolFromSmiles(p) for p in products_smiles.split(".") if p]

                # Check for silyl protection (alcohol to silyl ether)
                silyl_ether_patt = Chem.MolFromSmarts("[#6]-[#14](-[#6])(-[#6])-[#8]-[#6]")
                alcohol_patt = Chem.MolFromSmarts("[#8H]-[#6]")

                # Protection: alcohol in reactants, silyl ether in products
                if any(mol.HasSubstructMatch(alcohol_patt) for mol in reactants if mol) and any(
                    mol.HasSubstructMatch(silyl_ether_patt) for mol in products if mol
                ):
                    protection_found = True
                    print("Found silyl protection reaction")

                # Deprotection: silyl ether in reactants, alcohol in products
                if any(mol.HasSubstructMatch(silyl_ether_patt) for mol in reactants if mol) and any(
                    mol.HasSubstructMatch(alcohol_patt) for mol in products if mol
                ):
                    deprotection_found = True
                    print("Found silyl deprotection reaction")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found
