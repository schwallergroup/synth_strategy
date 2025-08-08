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
    This function detects amide bond formation in the synthetic route.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide bond formation
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6][#6](=[O])[O;H1]")
            amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O)]")
            amide_pattern = Chem.MolFromSmarts("[#6][#6](=[O])[#7]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            has_acid = any(
                mol and mol.HasSubstructMatch(carboxylic_acid_pattern)
                for mol in reactant_mols
                if mol
            )
            has_amine = any(
                mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
            )
            has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

            if has_acid and has_amine and has_amide:
                print("Found amide bond formation")
                amide_formation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
