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
    This function detects phenol alkylation with a haloalkane to form an aryl ether.
    """
    phenol_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for phenol and haloalkane in reactants
            phenol_pattern = Chem.MolFromSmarts("[OH]c")
            haloalkane_pattern = Chem.MolFromSmarts("[C][Cl,Br,I]")
            aryl_ether_pattern = Chem.MolFromSmarts("cO[C]")

            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and aryl_ether_pattern:
                has_phenol = any(
                    r and r.HasSubstructMatch(phenol_pattern) for r in reactants_mols if r
                )
                has_haloalkane = any(
                    r and r.HasSubstructMatch(haloalkane_pattern) for r in reactants_mols if r
                )
                product_has_aryl_ether = product_mol.HasSubstructMatch(aryl_ether_pattern)

                if has_phenol and has_haloalkane and product_has_aryl_ether:
                    phenol_alkylation_detected = True
                    print("Phenol alkylation with haloalkane detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return phenol_alkylation_detected
