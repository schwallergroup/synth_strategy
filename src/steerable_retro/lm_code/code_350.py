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
    This function detects if the synthetic route involves an etherification reaction
    between a phenol and an alkyl halide.
    """
    has_etherification = False

    def dfs_traverse(node):
        nonlocal has_etherification

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_pattern = "c[OH]"
                # Check for alkyl halide in reactants
                alkyl_halide_pattern = "[CX4][Cl,Br,I]"
                # Check for ether in product
                ether_pattern = "cO[CX4]"

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    has_phenol = reactants_mol.HasSubstructMatch(Chem.MolFromSmarts(phenol_pattern))
                    has_alkyl_halide = reactants_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(alkyl_halide_pattern)
                    )
                    has_ether = product_mol.HasSubstructMatch(Chem.MolFromSmarts(ether_pattern))

                    if has_phenol and has_alkyl_halide and has_ether:
                        print("Detected etherification between phenol and alkyl halide")
                        has_etherification = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_etherification
