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
    This function detects a synthetic strategy involving protection of a carboxylic acid
    with a tert-butyl group.
    """
    protection_found = False

    def dfs_traverse(node):
        nonlocal protection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
                # Check for tert-butyl ester in product
                tbutyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")

                has_acid = False

                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(carboxylic_acid_pattern):
                            has_acid = True
                            break
                    except:
                        continue

                try:
                    p_mol = Chem.MolFromSmiles(product)
                    if has_acid and p_mol and p_mol.HasSubstructMatch(tbutyl_ester_pattern):
                        print("Found carboxylic acid protection with tert-butyl group")
                        protection_found = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_found
