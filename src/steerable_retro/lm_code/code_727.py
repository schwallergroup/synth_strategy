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
    This function detects if the synthetic route involves alkylation of a phenol
    to form an alkyl aryl ether.
    """
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    alkyl_aryl_ether_pattern = Chem.MolFromSmarts("cO[CH]")

    phenol_alkylation = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check for phenol in reactants
                    reactants_have_phenol = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(phenol_pattern)
                        for r in reactants_smiles
                        if Chem.MolFromSmiles(r)
                    )

                    # Check for alkyl aryl ether in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if (
                        reactants_have_phenol
                        and product_mol
                        and product_mol.HasSubstructMatch(alkyl_aryl_ether_pattern)
                    ):
                        phenol_alkylation = True
                        print("Detected phenol alkylation to form alkyl aryl ether")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return phenol_alkylation
