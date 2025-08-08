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
    This function detects if the synthesis route involves reduction of nitro group to amine.
    """
    nitro_pattern = Chem.MolFromSmarts("[c]-[N+](=[O])-[O-]")
    aniline_pattern = Chem.MolFromSmarts("[c]-[NH2]")
    nitro_reduction = False

    def dfs_traverse(node):
        nonlocal nitro_reduction

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction: ArNO2 -> ArNH2
                reactant_has_nitro = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                        reactant_has_nitro = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                product_has_amine = product_mol and product_mol.HasSubstructMatch(aniline_pattern)

                if reactant_has_nitro and product_has_amine:
                    print("Nitro reduction to amine detected")
                    nitro_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction
