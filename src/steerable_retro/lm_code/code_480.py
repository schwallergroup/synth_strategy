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
    Detects if the synthesis involves coupling with a protected amine fragment
    that is later deprotected.
    """
    protected_amine_found = False
    subsequent_deprotection = False

    def dfs_traverse(node):
        nonlocal protected_amine_found, subsequent_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection
            boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])[#8]C(=O)[#7]")

            # Check if product has Boc group
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                print("Protected amine (Boc) detected in product")
                protected_amine_found = True

            # Check for deprotection (Boc group in reactant but not in product)
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and reactant_mol.HasSubstructMatch(boc_pattern)
                    and product_mol
                    and not product_mol.HasSubstructMatch(boc_pattern)
                ):
                    print("Boc deprotection detected")
                    subsequent_deprotection = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protected_amine_found and subsequent_deprotection
