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
    This function detects if the synthesis involves formation of an aryl ether.
    """
    aryl_ether_pattern = Chem.MolFromSmarts("c[#8][#6]")
    phenol_pattern = Chem.MolFromSmarts("c[#8H]")

    aryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal aryl_ether_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and all(r for r in reactant_mols):
                # Check if any reactant has a phenol and product has an aryl ether
                if any(
                    r.HasSubstructMatch(phenol_pattern) for r in reactant_mols if r
                ) and product_mol.HasSubstructMatch(aryl_ether_pattern):
                    aryl_ether_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Aryl ether formation strategy: {aryl_ether_formation}")
    return aryl_ether_formation
