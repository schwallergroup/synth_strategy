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
    This function detects if the synthetic route involves an amide bond disconnection.
    """
    found_disconnection = False

    def dfs_traverse(node):
        nonlocal found_disconnection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains amide bond that's broken in reactants
            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                amide_pattern = Chem.MolFromSmarts("[NH][C](=[O])")
                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if we have multiple fragments in reactants
                    if len(reactants) > 1:
                        # Check if one fragment has carboxylic acid and another has amine
                        has_acid = False
                        has_amine = False

                        for r in reactants:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol:
                                if r_mol.HasSubstructMatch(Chem.MolFromSmarts("[OH][C](=[O])")):
                                    has_acid = True
                                if r_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                                    has_amine = True

                        if has_acid and has_amine:
                            print("Found amide bond disconnection")
                            found_disconnection = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_disconnection
