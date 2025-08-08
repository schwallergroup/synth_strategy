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
    Detects the strategy of constructing a benzoxazinone scaffold followed by
    transformation to quinazolinone scaffold in early synthesis steps.
    """
    # Initialize tracking variables
    benzoxazinone_formed = False
    quinazolinone_formed = False
    benzoxazinone_to_quinazolinone = False

    # SMARTS patterns
    benzoxazinone_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1[#6](=[#8])[#8][#6]")
    quinazolinone_pattern = Chem.MolFromSmarts(
        "[#6]1[#7][#6]2[#6][#6][#6][#6][#6]2[#6](=[#8])[#7]1"
    )

    def dfs_traverse(node, depth):
        nonlocal benzoxazinone_formed, quinazolinone_formed, benzoxazinone_to_quinazolinone

        if node["type"] == "mol":
            # Check molecule for scaffolds
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if benzoxazinone_pattern and mol.HasSubstructMatch(benzoxazinone_pattern):
                    print(f"Benzoxazinone detected at depth {depth}")
                    benzoxazinone_formed = True

                if quinazolinone_pattern and mol.HasSubstructMatch(quinazolinone_pattern):
                    print(f"Quinazolinone detected at depth {depth}")
                    quinazolinone_formed = True

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if this reaction transforms benzoxazinone to quinazolinone
            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(
                    benzoxazinone_pattern
                ) and product_mol.HasSubstructMatch(quinazolinone_pattern):
                    print(
                        f"Benzoxazinone to quinazolinone transformation detected at depth {depth}"
                    )
                    benzoxazinone_to_quinazolinone = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route, 0)

    # Check if the strategy is present
    strategy_present = (
        benzoxazinone_formed and quinazolinone_formed and benzoxazinone_to_quinazolinone
    )
    print(f"Heterocyclic scaffold construction strategy detected: {strategy_present}")
    return strategy_present
