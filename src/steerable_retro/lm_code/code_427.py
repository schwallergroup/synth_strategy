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
    This function detects if the synthetic route involves an amine to azide transformation.
    """
    found_transformation = False

    def dfs_traverse(node):
        nonlocal found_transformation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains azide and reactant contains primary amine
            product_mol = Chem.MolFromSmiles(product)
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants]

            azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=N")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            if product_mol and any(Chem.MolFromSmiles(r) for r in reactants):
                has_azide = product_mol.HasSubstructMatch(azide_pattern)
                has_amine = any(r.HasSubstructMatch(amine_pattern) for r in reactants_mols if r)

                if has_azide and has_amine:
                    print("Found amine to azide transformation")
                    found_transformation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_transformation
