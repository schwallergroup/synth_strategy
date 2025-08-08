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
    This function detects if the synthesis involves amide bond formation.
    """
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                amine_pattern = Chem.MolFromSmarts("[NH]")
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if any reactant has carboxylic acid
                has_carboxylic_acid = any(
                    r is not None and r.HasSubstructMatch(carboxylic_acid_pattern)
                    for r in reactants
                )

                # Check if any reactant has amine
                has_amine = any(
                    r is not None and r.HasSubstructMatch(amine_pattern) for r in reactants
                )

                # Check if product has amide
                has_amide_in_product = product is not None and product.HasSubstructMatch(
                    amide_pattern
                )

                # Check if this is an amide formation reaction
                if has_carboxylic_acid and has_amine and has_amide_in_product:
                    print("Found amide bond formation")
                    has_amide_formation = True
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_amide_formation
