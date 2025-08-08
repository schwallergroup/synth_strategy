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
    This function detects a strategy involving the transformation of an amine to a thiazole
    through isothiocyanate and thiourea intermediates.
    """
    amine_present = False
    isothiocyanate_present = False
    thiourea_present = False
    thiazole_present = False

    def dfs_traverse(node, depth=0):
        nonlocal amine_present, isothiocyanate_present, thiourea_present, thiazole_present

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for functional groups in reactants and products
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for amine (typically at highest depth)
                    if depth >= 4:
                        amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
                        if product_mol.HasSubstructMatch(amine_pattern):
                            amine_present = True
                            print(f"Found amine at depth {depth}")

                    # Check for isothiocyanate
                    isothiocyanate_pattern = Chem.MolFromSmarts("[#6]-[N]=[C]=[S]")
                    if product_mol.HasSubstructMatch(isothiocyanate_pattern):
                        isothiocyanate_present = True
                        print(f"Found isothiocyanate at depth {depth}")

                    # Check for thiourea
                    thiourea_pattern = Chem.MolFromSmarts("[#6]-[NH]-[C](=[S])-[NH2]")
                    if product_mol.HasSubstructMatch(thiourea_pattern):
                        thiourea_present = True
                        print(f"Found thiourea at depth {depth}")

                    # Check for thiazole
                    thiazole_pattern = Chem.MolFromSmarts("c1nc([#6])cs1")
                    if product_mol.HasSubstructMatch(thiazole_pattern):
                        thiazole_present = True
                        print(f"Found thiazole at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amine_present and isothiocyanate_present and thiourea_present and thiazole_present
