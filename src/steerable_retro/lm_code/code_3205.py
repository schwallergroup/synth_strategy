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
    Detects a strategy involving amide formation from carboxylic acid and amine.
    """
    found_amide_disconnection = False

    def dfs_traverse(node):
        nonlocal found_amide_disconnection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[C](=[O])-[NH]-[C]")
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])-[OH]")
                amine_pattern = Chem.MolFromSmarts("[NH2]-[C]")

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check if product contains amide bond
                    if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                        # Check if reactants include carboxylic acid and amine
                        has_acid = False
                        has_amine = False

                        for r in reactants:
                            try:
                                r_mol = Chem.MolFromSmiles(r)
                                if r_mol and r_mol.HasSubstructMatch(carboxylic_acid_pattern):
                                    has_acid = True
                                if r_mol and r_mol.HasSubstructMatch(amine_pattern):
                                    has_amine = True
                            except:
                                continue

                        if has_acid and has_amine:
                            found_amide_disconnection = True
                            print(f"Found amide formation from acid and amine: {rsmi}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_amide_disconnection
