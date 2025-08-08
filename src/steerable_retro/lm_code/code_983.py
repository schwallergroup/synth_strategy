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
    This function detects if the synthetic route involves esterification
    of a carboxylic acid to form a methyl ester.
    """
    esterification_detected = False

    def dfs_traverse(node):
        nonlocal esterification_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
            methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check if product has methyl ester
                if product_mol and product_mol.HasSubstructMatch(methyl_ester_pattern):
                    # Check if any reactant has carboxylic acid
                    for reactant_smiles in reactants_smiles:
                        try:
                            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                carboxylic_acid_pattern
                            ):
                                esterification_detected = True
                                print("Detected esterification of carboxylic acid to methyl ester")
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
    return esterification_detected
