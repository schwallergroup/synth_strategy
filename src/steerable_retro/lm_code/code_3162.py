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
    This function detects if a BOC protection/deprotection strategy is used
    for a piperazine nitrogen.
    """
    boc_protected_piperazine = False
    boc_deprotection = False

    def dfs_traverse(node):
        nonlocal boc_protected_piperazine, boc_deprotection

        if node["type"] == "mol":
            if "smiles" in node:
                smiles = node["smiles"]

                # Check for BOC-protected piperazine
                boc_piperazine_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]1[C][C][N][C][C]1")
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol and mol.HasSubstructMatch(boc_piperazine_pattern):
                        boc_protected_piperazine = True
                        print(f"Found BOC-protected piperazine: {smiles}")
                except:
                    pass

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for BOC deprotection
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")

                try:
                    # Check if any reactant has BOC group
                    reactant_has_boc = False
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(boc_pattern):
                            reactant_has_boc = True
                            break

                    # Check if product doesn't have BOC group
                    product_mol = Chem.MolFromSmiles(product)
                    product_has_boc = product_mol and product_mol.HasSubstructMatch(boc_pattern)

                    if reactant_has_boc and not product_has_boc:
                        boc_deprotection = True
                        print(f"Detected BOC deprotection: {rsmi}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if both BOC protection and deprotection are detected
    return boc_protected_piperazine and boc_deprotection
