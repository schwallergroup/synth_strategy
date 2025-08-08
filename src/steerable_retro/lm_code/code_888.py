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
    This function detects Wittig olefination.
    Looks for C=C bond formation using phosphorane and aldehyde.
    """
    wittig_detected = False

    def dfs_traverse(node):
        nonlocal wittig_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")

                # Check for phosphorane in reactants
                phosphorane_pattern = Chem.MolFromSmarts("P(c1ccccc1)(c1ccccc1)=C")

                # Check for C=C in product
                alkene_pattern = Chem.MolFromSmarts("C=C-C(=O)O")

                has_aldehyde = False
                has_phosphorane = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(aldehyde_pattern):
                            has_aldehyde = True
                        if mol and mol.HasSubstructMatch(phosphorane_pattern) or "P" in reactant:
                            has_phosphorane = True
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    has_alkene = product_mol and product_mol.HasSubstructMatch(alkene_pattern)

                    if has_aldehyde and has_phosphorane and has_alkene:
                        print("Detected Wittig olefination")
                        wittig_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return wittig_detected
