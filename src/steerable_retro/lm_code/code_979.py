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
    This function detects if the synthetic route involves benzimidazole formation
    from an ortho-nitroaniline precursor.
    """
    benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2[nH]1")
    benzimidazole_formed = False

    def dfs_traverse(node):
        nonlocal benzimidazole_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants contain ortho-nitroaniline pattern
            ortho_nitroaniline_pattern = Chem.MolFromSmarts("[NH2]c1c([N+](=[O])[O-])cccc1")

            # Check if product contains benzimidazole
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(benzimidazole_pattern):
                    # Check if any reactant has ortho-nitroaniline pattern
                    for reactant_smiles in reactants_smiles:
                        try:
                            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                ortho_nitroaniline_pattern
                            ):
                                benzimidazole_formed = True
                                print("Detected benzimidazole formation from ortho-nitroaniline")
                                break
                        except:
                            continue
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return benzimidazole_formed
