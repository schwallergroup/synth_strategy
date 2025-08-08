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
    Detects if the synthetic route involves a Boc protection-deprotection sequence.
    """
    boc_protected_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_protected_found, deprotection_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc group in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
                    if product_mol.HasSubstructMatch(boc_pattern):
                        boc_protected_found = True

                # Check for Boc deprotection (Boc in reactants but not in product)
                reactant_has_boc = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
                        if reactant_mol.HasSubstructMatch(boc_pattern):
                            reactant_has_boc = True

                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
                    if reactant_has_boc and not product_mol.HasSubstructMatch(boc_pattern):
                        deprotection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(
        f"Boc protection-deprotection sequence detected: {boc_protected_found and deprotection_found}"
    )
    return boc_protected_found and deprotection_found
