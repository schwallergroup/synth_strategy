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
    This function detects if the synthetic route involves amide bond formation from a carboxylic acid and an amine.
    """
    amide_formation = False

    def dfs_traverse(node):
        nonlocal amide_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
            has_carboxylic_acid = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                    has_carboxylic_acid = True
                    break

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[NH2][#6]")
            has_amine = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                    has_amine = True
                    break

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH][#6]")
            product_mol = Chem.MolFromSmiles(product)

            if (
                has_carboxylic_acid
                and has_amine
                and product_mol
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                amide_formation = True
                print(f"Amide bond formation detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Amide bond formation strategy: {amide_formation}")
    return amide_formation
